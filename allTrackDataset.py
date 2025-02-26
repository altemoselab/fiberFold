# Classes to load data
# using pytorch DataLoader and DataSet

##Danilo Dubocanin ##

import sys 
import numpy as np
import pandas as pd
import pyBigWig as pbw
import os
from tqdm import tqdm
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from shapeData import GenomicFeatureSingleThread
from shapeData import HiCFeature
import pyfaidx
import pyfastx 
import gzip



def main():
    
    data = GenomeDataset(sys.argv[1])

    dataloader = DataLoader(data,shuffle=True,batch_size=8)


class GenomeDataset(Dataset):
    '''
    Load all chromosomes
    '''
    def __init__(self, datasheet, ma_bw, mc_bw, ctcf,ctcf_dir, hic=None ,rev=True, shift = False, noise=True):
        with gzip.open(datasheet, 'rt') as handle:
            self.files = handle.read().splitlines()
        self.rev = rev
        self.ma_bw_dir = ma_bw
        self.mc_bw_dir = mc_bw
        self.ctcf = ctcf
        self.ctcf_dir = ctcf_dir
        self.hic_dir = hic
        self.shift = shift
        self.noise = noise

    def __getitem__(self, idx):

        coords = str(self.files[idx])

        chrom, start, stop  = coords.split("\t")
        start, stop  = int(start), int(stop)
        if self.shift == True:
            
            max_shift = int(0.005 * (stop-start))
            shift_val = int( np.random.randint(-max_shift,max_shift) )

            start = int(start + shift_val)
            stop = int(stop + shift_val)

        self.ma_bw = GenomicFeatureSingleThread(f'{self.ma_bw_dir}/{chrom}.bw',norm=False) 
        self.mc_bw = GenomicFeatureSingleThread(f'{self.mc_bw_dir}/{chrom}.bw',norm=False)  
        self.ctcf_bw = GenomicFeatureSingleThread(f'{self.ctcf}/{chrom}.bw',norm=False) # for actual processing
        self.ctcf_dir_bw = GenomicFeatureSingleThread(f'{self.ctcf_dir}/{chrom}.bw',norm=False) # for actual processing
        # self.ctcf = GenomicFeatureSingleThread('/Users/danilodubocanin/work/altemose/fiber2ctcf/test_data/ctcf_preds_all_v3.bw',norm=False)



        # calculate mA


        mA_feature = self.ma_bw.get(chrom,start,stop,add_noise=self.noise)

        #calculate mC 
        mC_feature = self.mc_bw.get(chrom,start,stop,add_noise=self.noise)
        # calculate CTCF Chip with log scale
        ctcf_feature = self.ctcf_bw.get(chrom,start,stop,add_noise=self.noise)
        ctcf_dir_feature = self.ctcf_dir_bw.get(chrom,start,stop,add_noise=self.noise)


        # inputs
        input_features = np.vstack([mA_feature, mC_feature, ctcf_feature, ctcf_dir_feature]).T


        # output feature will be Hi-C map here
        if self.hic_dir:
            hic_map = HiCFeature(f'{self.hic_dir}/{chrom}.npz').get(start)

            if self.rev == True:

                input_features, matrix = self.flip_seq(input_features,hic_map)

            else:
                matrix = hic_map

            return input_features.astype(np.float32), matrix.astype(np.float32), tuple([chrom,start,stop])

        else:

            return input_features.astype(np.float32), tuple([chrom,start,stop])


    def flip_seq(self, features, mat):


        if np.random.rand(1) < 0.5 : 


            features_r = np.fliplr(features)

            mat_r = np.flip(mat)

        else:

            features_r, mat_r =  features, mat


        return features_r, mat_r
    
    def __len__(self):
        return len(self.files)
        
    
if __name__ == '__main__':
    main()
