import argparse 
import sys
import torch
import numpy as np
import model_simple as models
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from allTrackDataset import GenomeDataset
from tqdm import tqdm
import pandas as pd
import pickle 


def prediction(output_path, model_path, datasheet, ma,mc,ctcf_fp,ctcf_dir):
	
	### initialize model 	####

	# define model name 
	model_name = 'ConvTransModel'

	# get attributes of that model type
	ModelClass = getattr(models,model_name)

	# initialize model with num_genomic_features
	model = ModelClass(4)

	# load the trained model
	model = load_checkpoint(model,model_path)

	all_predictions = []
	device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
	# iterate through the data
	with torch.no_grad():
		for d in tqdm(getData(datasheet, ma, mc, ctcf_fp, ctcf_dir)):

			inputs,chrom_info = d
			inputs = inputs.to(device)
			pred = model(inputs).detach().cpu().numpy()
			
			all_predictions.append((pred, chrom_info))
			torch.cuda.empty_cache()


	with open(output_path, 'wb') as handle:
		pickle.dump(all_predictions, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_checkpoint(model, model_path):
    print('Loading weights')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    checkpoint = torch.load(model_path, map_location=device)
    model_weights = checkpoint['state_dict']

    # Edit keys
    for key in list(model_weights):
        model_weights[key.replace('model.', '')] = model_weights.pop(key)
    model.load_state_dict(model_weights)
    model.eval()

    return model

def getData(datasheet, ma, mc, ctcf,ctcf_dir):
	dl = DataLoader(GenomeDataset(datasheet, ma, mc, ctcf,ctcf_dir, rev = False),batch_size=1,
																		num_workers=4,
																		pin_memory=True,
																		prefetch_factor=1,
																		persistent_workers=True
																		)
	for step, t in enumerate(dl):
		# continue
		yield proc_batch(t)

def proc_batch(batch):
    
	features , chrom_info = batch

	return features, chrom_info

def main():

	parser = argparse.ArgumentParser(description='fiber2hic 2M prediction module')

	parser.add_argument('--out', dest='output_path', default='outputs',
				         help='output path for storing results')

	parser.add_argument('--model', dest='model_path', 
	                    help='Path to the model checkpoint to use for predictions', required=True)
	parser.add_argument('--ds', dest='datasheet',
							help = 'Path To datasheet with files')
	
	parser.add_argument('--mA',dest='ma_bw', default='data',
	                    help='mA bigWig', required=True)
	parser.add_argument('--mC',dest='mc_bw', default='data',
	                    help='mC_bigWig', required=True)
	parser.add_argument('--ctcf_dir',dest='ctcf_dir', default='data',
	                    help='ctcf direction bigwig', required=True)
	parser.add_argument('--ctcf_fp',dest='ctcf_fp', default='data',
	                    help='ctcf fiberseq footprints scores', required=True)

	args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])


	prediction(args.output_path,args.model_path, args.datasheet, args.ma_bw, args.mc_bw, args.ctcf_fp, args.ctcf_dir)

if __name__ == '__main__':
    main()