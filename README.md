FiberFold enables the prediction of Hi-C maps (2MB window, 10kb resolution) from a single assay, Fiber-seq.

First, generate bigWigs of FIRE, mCpG, CTCF footprinting, and CTCF motif direction.

Second, split these data by chromsosome and store all chroms from a single feature track in a directory.

Third, generate your prediction windows. Must be 2097152 bases long interval, organized as a bed file.

Fourth, predict with predict.py or alternatively, train a model with main.py. 

Example directory structure should look like (in file_prefix.bw, file prefix must be chromosome name):
data/
  cpg/
    chr1.bw
    chr2.bw
  fire/
    chr1.bw
    chr2.bw
  ctcf_footprint/
    chr1.bw
    chr2.bw
  ctcf_direction/
    chr1.bw
    chr2.bw


To generate phased Hi-C maps, phase Fiber-seq reads prior to the first step.

This code-base and model architecture is inspired by, and largely a fork of, [C.Origami](https://github.com/tanjimin/C.Origami)  (Tan et. al. *Nature Biotechnology* 2023, https://www.nature.com/articles/s41587-022-01612-8).

To generate FIRES and CTCF footprinting scores, please use [fiberTools](https://github.com/fiberseq/fibertools-rs) (Jha & Bohaczuk et al. *Genome Research*, 2024)


