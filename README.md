FiberFold enables the prediction of Hi-C maps (2MB window, 10kb resolution) from a single assay, Fiber-seq.

First, generate bigWigs of FIRE, mCpG, CTCF footprinting, and CTCF motif direction.
Second, split these data by chromsosome and store all chroms from a single feature track in a directory. 
Third, generate your prediction windows. Must be 2097152 bases long interval, organized as a bed file.
Fourth, predict with predict.py or alternatively, train a model with main.py. 

To generate phased Hi-C maps, phase Fiber-seq reads prior to the first step.

This code-base and model architecture is inspired by, and largely a fork of, C.Origami (Tan et. al. Nature Biotech 2023).


