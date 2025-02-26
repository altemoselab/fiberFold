import sys 
import numpy as np
# requires sorted bed input 


outputs = []

with open(sys.argv[1],'r') as handle:

	last_out = None

	for line in handle:

		chrom,start,stop,name,score,strand = line.strip().split()

		if last_out is None:

			outputs.append(line.strip())

			last_out = ((chrom,start,stop,name,score,strand ))

		else:

			if int(start) - int(last_out[2]) <= 0:
				score_idx = np.argmax([float(last_out[-2]) , float(score)])

				if score_idx == 1:
					outputs[-1] = line.strip()
					last_out = ((chrom,start,stop,name,score,strand ))
				else:
					continue
			else:
				outputs.append(line.strip())

				last_out = ((chrom,start,stop,name,score,strand ))


for x in outputs:
	print(x)