## FiberFold enables the prediction of Hi-C maps (2MB window, 10kb resolution) from a single assay, Fiber-seq (PacBio).

First, generate bigWigs of FIRE, mCpG, CTCF footprinting, and CTCF motif direction.
  (pb-CpG-tools)[https://github.com/PacificBiosciences/pb-CpG-tools] for generating CpG (default on PacBio output)
  [fiberTools](https://github.com/fiberseq/fibertools-rs) FIRE for generating accessible regulatory element track
  [fiberTools](https://github.com/fiberseq/fibertools-rs) footprint for generating CTCF footprinting information, followed by misc/encodeFP.py + misc/filterOverlapMotif.py. CTCF Motif to use: MA0139.1.
  misc/encodeMotifDirection.py to encode CTCF motif directionality.
 
Second, split these data by chromsosome and store all chroms from a single feature track in a directory.
  Example directory structure should look like (in file_prefix.bw, file prefix must be chromosome name):
  
      data/
        cpg/
          chr1.bw
          chr2.bw
          ...
        fire/
          chr1.bw
          chr2.bw
          ...
        ctcf_footprint/
          chr1.bw
          chr2.bw
          ...
        ctcf_motif_direction/
          chr1.bw
          chr2.bw
          ...
Third, generate your prediction windows. Must be 2097152 bases long interval, organized as a bed file.

Fourth, predict with predict.py or alternatively, train a model with main.py. 

To generate phased Hi-C maps, phase Fiber-seq reads prior to the first step.

Attributions:
The model architecture is inspired by, and largely a fork of, [C.Origami](https://github.com/tanjimin/C.Origami)  (Tan et. al. *Nature Biotechnology* 2023, https://www.nature.com/articles/s41587-022-01612-8).


