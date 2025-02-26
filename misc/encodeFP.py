import numpy as np
from tqdm import tqdm
import pyBigWig
import gzip
from pyfaidx import Fasta

# to Run:
# python3 encodeFP.py input_fp.gz chromsizes genome.fa output.bw 



fp_file = sys.argv[1]

with open(sys.argv[2],'r') as handle:
    chromsizes = [(c.strip().split()[0],int(c.strip().split()[1])) for c in handle.readlines() if "_" not in c.split()[0]]

# need genome to quickly lookup specific CTCF motif sequence
genome = Fasta(sys.argv[3],sequence_always_upper=True,one_based_attributes=False)

out_file = sys.argv[4]



def processMotif(line,):
    
    split_line = line.strip().split()
    if int(split_line[4]) == 0:
        return None,None
    modules = [[0,6],[6,13],[13,19]]
    
    chrom, start, end, strand, n_spanning_fibers, n_spanning_msps, n_overlapping_nucs, m1, m2, m3, fp_codes, fire_quals, fiber_names = line.strip().split()
    
    seq = genome[chrom][int(start):int(end)].seq
    
    
    relevant_modules = []

    for c,m in enumerate(modules):
        if "A" in seq[m[0]:m[1]] or "T" in seq[m[0]:m[1]]:
            relevant_modules.append(c + 1)
    # first need to calculate the number of motifs footprinted 
    # check each motif individually to see which modules can even be footprinted (if there is A/T )
    
    true_footprints = 0
    
    for fp_code in fp_codes.split(','):
        
        is_fp = processCode(int(fp_code), relevant_modules)
        if is_fp == True: true_footprints += 1
    
    return true_footprints/int(n_spanning_fibers), n_spanning_fibers
    
    
def processCode(fp_code,relevant_modules):
    
    coded = []
    if (fp_code & 1) > 0:
        coded.append(True)
    else:
        coded.append(False)
    
    for m in relevant_modules:
        
        if (fp_code & (1 << m)) > 0:
            coded.append(True)
        else:
            coded.append(False)
            
    if len(coded) == np.sum(coded):
        return True
    else:
        return False



bw = pyBigWig.open(out_file,'w')
bw.addHeader(chromsizes)


for c in chromsizes[:-3]:
    chrom=c[0]
    with gzip.open(fp_file,'rt') as handle:
        all_out_vals = []
        all_out_starts = []
        all_out_chroms = []
    #     current_chrom = None
        for line in tqdm(handle):

            if 'chrY' in line:
                continue
            if 'chrM' in line:
                continue
            if 'chrEBV' in line:
                continue
            if "#" in line: 
                continue
               
            split_line = line.split()
            if split_line[0] != chrom: continue
            motif_vals, cov = processMotif(line)

            if type(motif_vals) == type(None):
                continue

            all_out_vals.append(float(motif_vals))
            all_out_starts.append(int(split_line[1]))

            all_out_chroms.append(str(split_line[0]))

        all_out_ends = [int(s+19) for s in all_out_starts]

        ns = np.array(all_out_starts)
        ne = np.array(all_out_ends)
        nv = np.array(all_out_vals)
        nc = np.array(all_out_chroms)
        bw.addEntries(list(nc),list(ns.astype(int)),list(ne.astype(int)),list(nv.astype(float)))
            
            
bw.close()        



