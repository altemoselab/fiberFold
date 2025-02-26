import numpy as np 
import pickle as pkl 
import matplotlib.pyplot as plt
import pyBigWig as bw 
from tqdm import tqdm
from pybedtools import BedTool


# usage: python3 encodeMotifDirection.py motif_file.bed chromsizes out.bw

motif_file = BedTool(sys.argv[1])

with open(sys.argv[2],'r') as handle:
    chromsizes = [(c.strip().split()[0],int(c.strip().split()[1])) for c in handle.readlines() if "_" not in c.split()[0]]

out_file = bw.open(sys.argv[3],'w')
out_file.addHeader(chromsizes)


def rle(inarray): 
	""" run length encoding. Partial credit to R rle function. 
	Multi datatype arrays catered for including non Numpy
	returns: tuple (runlengths, startpositions, values) 
	found here: https://stackoverflow.com/questions/1066758/
	find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi"""
	ia = np.asarray(inarray)                # force numpy
	n = len(ia)

	y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
	i = np.append(np.where(y), n - 1)   # must include last element posi
	z = np.diff(np.append(-1, i))       # run lengths
	p = np.cumsum(np.append(0, z))[:-1] # positions
	return(z, p, ia[i])


for chrom,chrom_size in tqdm(chromsizes):

	chrom_bed = BedTool(f'{chrom}\t1\t{chrom_size}', from_string=True)
	    
	val_array   = np.zeros(shape=int(chrom_size),   dtype = np.int8 )

	for motif in motif_file.intersect(chrom_bed,wa=True,u=True):

		motif_start, motif_stop = int(motif.start), int(motif.stop)

		if motif[-1] == "+":
			val = 1
		elif motif[-1] == "-":
			val = -1 
		else:
			val = 0

		val_array[motif_start:motif_stop] = val

	rle_arr = rle(val_array)
	rle_arr = np.vstack(rle_arr)

	c = list([str(chrom)]*(len(rle_arr[0])))
	s = list(rle_arr[1].astype(int))       
	e = list(np.add(rle_arr[1],rle_arr[0]).astype(int))
	v = list(rle_arr[2].astype(float))

	out_file.addEntries(c,s,ends=e,values=v)

out_file.close()





