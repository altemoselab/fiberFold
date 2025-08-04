import sys
import torch
import argparse
import numpy as np
import pytorch_lightning as pl
import pytorch_lightning.callbacks as callbacks
from pl_bolts.optimizers.lr_scheduler import LinearWarmupCosineAnnealingLR

import model_simple
from allTrackDataset import GenomeDataset


def main():
    args = init_parser()
    init_training(args)

def init_parser():
  parser = argparse.ArgumentParser(description='C.Origami Training Module.')

  # Data and Run Directories
  parser.add_argument('--seed', dest='run_seed', default=2806,
                        type=int,
                        help='Random seed for training')
  parser.add_argument('--save_path', dest='run_save_path', default='checkpoints',
                        help='Path to the model checkpoint')

  # Data directories
  parser.add_argument('--datasheet_val', dest='datasheet_val', default='data',
                        help='datasheet with al data locations', required=True)

  parser.add_argument('--datasheet_train', dest='datasheet_train', default='data',
                        help='datasheet with al data locations', required=True)

  parser.add_argument('--mA',dest='ma_bw', default='data',
                        help='mA bigWig', required=True)
  parser.add_argument('--mC',dest='mc_bw', default='data',
                        help='mC_bigWig', required=True)
  parser.add_argument('--hic',dest='hic', default='data',
                        help='hic dir, npz format', required=True)
  parser.add_argument('--ctcf',dest='ctcf', default='data',
                        help='ctcf chip bw (signal p-value)', required=True)

  parser.add_argument('--ctcf_dir',dest='ctcf_dir', default=None,
                        help='ctcf direction', required=True)


  parser.add_argument('--ckpt',dest='ckpt', default=None,
                        help='ctcf chip bw (signal p-value)', required=False)


  parser.add_argument('--n_feat',dest='n_feat', default=3,
                        help='number of features', required=False)
  # Train single stranded or both 

  args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
  return args

def init_training(args):

    # Early_stopping
    early_stop_callback = callbacks.EarlyStopping(monitor='val_loss', 
                                        min_delta=0.00, 
                                        patience=15,
                                        verbose=False,
                                        mode="min")
    # Checkpoints
    checkpoint_callback = callbacks.ModelCheckpoint(dirpath=f'{args.run_save_path}/models',
                                        save_top_k=10, 
                                        monitor='val_loss')

    # LR monitor
    lr_monitor = callbacks.LearningRateMonitor(logging_interval='epoch')

    # Logger
    csv_logger = pl.loggers.CSVLogger(save_dir = f'{args.run_save_path}/csv')
    all_loggers = csv_logger
    
    # Assign seed
    pl.seed_everything(args.run_seed, workers=True)
    pl_module = GANModule(args)
    pl_trainer = pl.Trainer(accelerator='gpu',strategy='ddp',
                            logger = all_loggers,
                            callbacks = [early_stop_callback,
                                         checkpoint_callback,
                                         lr_monitor],
                            max_epochs = 120, fast_dev_run=True) # fast_dev_run to see if model will initialize, remove when running entire model

    trainloader = pl_module.get_dataloader(args, 'train')
    valloader = pl_module.get_dataloader(args, 'val')

    pl_trainer.fit(pl_module, trainloader, valloader, ckpt_path = args.ckpt)


class TrainModule(pl.LightningModule):
    
    def __init__(self, args):
        super().__init__()
        self.model = self.get_model(args)
        self.args = args
        self.save_hyperparameters()

    def forward(self, x):
        return self.model(x)

    def proc_batch(self, batch):
        
        features, mat, chrom_info = batch

        return features, mat
    
    def training_step(self, batch, batch_idx):

        inputs, mat = self.proc_batch(batch)
        outputs = self(inputs)
        criterion = torch.nn.MSELoss()
        loss = criterion(outputs, mat)

        metrics = {'train_step_loss': loss}
        self.log_dict(metrics, batch_size = inputs.shape[0], prog_bar=True)

        return loss

    def validation_step(self, batch, batch_idx):
        ret_metrics = self._shared_eval_step(batch, batch_idx)
        return ret_metrics

    def test_step(self, batch, batch_idx):
        ret_metrics = self._shared_eval_step(batch, batch_idx)
        return ret_metrics

    def _shared_eval_step(self, batch, batch_idx):
            
        #original eval step

        inputs, mat = self.proc_batch(batch)
        outputs = self(inputs)

        criterion = torch.nn.MSELoss()
        #criterion = torch.nn.PoissonNLLLoss()
        loss = criterion(outputs, mat)

        return loss

    # Collect epoch statistics
    def training_epoch_end(self, step_outputs):
        step_outputs = [out['loss'] for out in step_outputs]
        ret_metrics = self._shared_epoch_end(step_outputs)
        metrics = {'train_loss' : ret_metrics['loss']
                  }
        self.log_dict(metrics, prog_bar=True)

    def validation_epoch_end(self, step_outputs):
        ret_metrics = self._shared_epoch_end(step_outputs)
        metrics = {'val_loss' : ret_metrics['loss']
                  }
        self.log_dict(metrics, prog_bar=True)

    def _shared_epoch_end(self, step_outputs):
        loss = torch.tensor(step_outputs).mean()
        return {'loss' : loss}

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), 
                                     lr = 4e-4,
                                     weight_decay = 0)

        scheduler = pl_bolts.optimizers.lr_scheduler.LinearWarmupCosineAnnealingLR(optimizer, 
            warmup_epochs=10, max_epochs=40)
        scheduler_config = {
            'scheduler': scheduler,
            'interval': 'epoch',
            'frequency': 1,
            'monitor': 'val_loss',
            'strict': True,
            'name': 'WarmupCosineAnnealing',
        }
        return {'optimizer' : optimizer, 'lr_scheduler' : scheduler_config}

    def get_dataset(self, args, mode):

        if mode == 'val':
            dataset = GenomeDataset(self.args.datasheet_val, self.args.ma_bw, self.args.mc_bw, self.args.ctcf, self.args.ctcf_dir, self.args.hic)
        elif mode == "train":
            dataset = GenomeDataset(self.args.datasheet_train, self.args.ma_bw, self.args.mc_bw, self.args.ctcf, self.args.ctcf_dir, self.args.hic,shift=True)
        
        return dataset

    def get_dataloader(self, args, mode):
        dataset = self.get_dataset(args, mode)

        if mode == 'train':
            shuffle = True
        else: # validation and test settings
            shuffle = False
        
        #### SET BATCH SIZE
        batch_size = 8
        num_workers = 8


        dataloader = torch.utils.data.DataLoader(
            dataset,
            shuffle=shuffle,
            batch_size=batch_size,

            num_workers=num_workers,
            pin_memory=True,
            prefetch_factor=1,
            persistent_workers=True
        )
        return dataloader

    def get_model(self, args):
        model_name = 'ConvTransModel' 
        num_genomic_features = int(args.n_feat)
        ModelClass = getattr(model_simple, model_name)
        model = ModelClass(num_genomic_features, mid_hidden = 256)
        return model


# TRAIN GAN 
class GANModule(pl.LightningModule):
    def __init__(self, args, lambda_adv=0.9):
        super().__init__()
        self.generator = self.get_model(args)
        if args.ckpt is not None:
            ckpt = torch.load(args.ckpt, map_location=self.device)
            if 'state_dict' in ckpt:
                pretrained_weights = {k.replace("model.", ""): v for k, v in ckpt['state_dict'].items() if "model." in k}
                self.generator.load_state_dict(pretrained_weights, strict=False)
            else:
                self.generator.load_state_dict(ckpt, strict=False)
            print(f"Loaded pretrained generator from: {args.ckpt}")


        self.discriminator = model_simple.Discriminator()
        self.args = args
        self.lambda_adv = lambda_adv
        self.save_hyperparameters()
        self.automatic_optimization = False

    def forward(self, x):
        return self.generator(x)

    def proc_batch(self, batch):
        features, mat, _ = batch
        return features, mat

    def training_step(self, batch, batch_idx):
        opt_g, opt_d = self.optimizers()
        inputs, real = self.proc_batch(batch)

        if self.current_epoch < 10:
            # Warm-up: only update discriminator
            fake = self.generator(inputs)
            opt_d.zero_grad()
            real_disc = self.discriminator(real.unsqueeze(1))
            fake_disc = self.discriminator(fake.detach().unsqueeze(1))
            labels = torch.cat([torch.ones_like(real_disc), torch.zeros_like(fake_disc)], dim=0)
            preds = torch.cat([real_disc, fake_disc], dim=0)
            disc_loss = torch.nn.binary_cross_entropy_with_logits(preds, labels)
            self.manual_backward(disc_loss)
            opt_d.step()
            self.log('train_disc_loss', disc_loss, prog_bar=True)
            


        else:
            fake = self.generator(inputs)
            B, H, W = real.shape

            # ----------------------
            #  Generator step
            # ----------------------
            opt_g.zero_grad()

            mse_loss = torch.nn.mse_loss(fake, real)
            disc_out = self.discriminator(fake.unsqueeze(1))  # (B,1,H,W) → logits
            adv_loss = torch.nn.binary_cross_entropy_with_logits(disc_out, torch.ones_like(disc_out))
            gen_loss = (self.lambda_adv * mse_loss) + ((1 - self.lambda_adv) * adv_loss)

            self.manual_backward(gen_loss)
            torch.nn.utils.clip_grad_norm_(self.generator.parameters(), max_norm=1.0)
            opt_g.step()

            # ----------------------
            #  Discriminator step
            # ----------------------
            opt_d.zero_grad()

            real_disc = self.discriminator(real.unsqueeze(1))
            fake_disc = self.discriminator(fake.detach().unsqueeze(1))

            labels = torch.cat([torch.ones_like(real_disc), torch.zeros_like(fake_disc)], dim=0)
            preds = torch.cat([real_disc, fake_disc], dim=0)

            disc_loss = self.discriminator.loss(preds, labels)

            self.manual_backward(disc_loss)
            opt_d.step()

            self.log_dict({
                'train_gen_loss': gen_loss,
                'train_disc_loss': disc_loss,
                'train_mse_loss': mse_loss,
                'train_adv_loss': adv_loss,
            }, prog_bar=True)

    def validation_step(self, batch, batch_idx):
        inputs, mat = self.proc_batch(batch)
        fake = self(inputs)
        mse_loss = torch.nn.functional.mse_loss(fake, mat)

        # Adversarial signal: is the discriminator fooled?
        fake_score = self.discriminator(fake.unsqueeze(1)).sigmoid().mean()

        self.log_dict({
            'val_mse_loss': mse_loss,
            'val_disc_conf_fake': fake_score,
        }, prog_bar=True)
        return mse_loss
    def configure_optimizers(self):
        opt_g = torch.optim.Adam(self.generator.parameters(), lr=4e-4)
        opt_d = torch.optim.Adam(self.discriminator.parameters(), lr=4e-4)

        scheduler_g = {
            "scheduler": LinearWarmupCosineAnnealingLR(opt_g, warmup_epochs=10, max_epochs=120),
            "interval": "epoch",
            "name": "gen_lr"
        }
        scheduler_d = {
            "scheduler": LinearWarmupCosineAnnealingLR(opt_d, warmup_epochs=10, max_epochs=120),
            "interval": "epoch",
            "name": "disc_lr"
        }

        return [opt_g, opt_d], [scheduler_g, scheduler_d]
    
    def get_model(self, args):
        model_name = 'ConvTransModel' 
        num_genomic_features = int(args.n_feat)
        ModelClass = getattr(model_simple, model_name)
        model = ModelClass(num_genomic_features, mid_hidden = 256)
        return model


    def get_dataset(self, args, mode):
        if mode == 'val':
            return GenomeDataset(args.datasheet_val, args.ma_bw, args.mc_bw, args.ctcf, args.ctcf_dir, args.hic)
        else:
            return GenomeDataset(args.datasheet_train, args.ma_bw, args.mc_bw, args.ctcf, args.ctcf_dir, args.hic, shift=True)

    def get_dataloader(self, args, mode):
        dataset = self.get_dataset(args, mode)

        if mode == 'train':
            shuffle = True
        else: # validation and test settings
            shuffle = False
        
        #### SET BATCH SIZE
        batch_size = 8
        num_workers = 8


        dataloader = torch.utils.data.DataLoader(
            dataset,
            shuffle=shuffle,
            batch_size=batch_size,

            num_workers=num_workers,
            pin_memory=True,
            prefetch_factor=1,
            persistent_workers=True
        )
        return dataloader

if __name__ == '__main__':
    main()

