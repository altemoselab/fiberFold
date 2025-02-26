## FiberFold enables the prediction of Hi-C maps (2MB window, 10kb resolution) from a single assay, Fiber-seq (PacBio).  
  
### First, generate bigWigs of FIRE, mCpG, CTCF footprinting, and CTCF motif direction.  
  [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) for generating CpG (default on PacBio output)  
  [fiberTools](https://github.com/fiberseq/fibertools-rs) FIRE for generating accessible regulatory element track  
  [fiberTools](https://github.com/fiberseq/fibertools-rs) footprint for generating CTCF footprinting information, followed by misc/encodeFP.py + misc/filterOverlapMotif.py. CTCF Motif to use: MA0139.1.  
  misc/encodeMotifDirection.py to encode CTCF motif directionality.   
  
  All tracks should be divided by depth to be in range 0-1.
   
### Second, split these data by chromsosome and store all chroms from a single feature track in a directory.  
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
            
### Third, generate your prediction target windows.   
  Must be 2097152 bases long interval, organized as a bed file. Can use [bedtools](https://github.com/arq5x/bedtools2) makewindows to do this across a genomic interval.
  
### Fourth, predict with predict.py or alternatively, train a model with main.py.   
  
To generate phased Hi-C maps, phase Fiber-seq reads prior to the first step.  

  To train:
  
      python main.py 
        --datasheet_val <alidation_chrom_windows.bed>
        --datasheet_train <training_chroms_windows.bed>
        --mA <fire/>
        --mC <cpg/>
        --ctcf_fp <ctcf_footprint>
        --ctcf_dir <ctcf_motif_direction>
        --n_feat 4
        --model model.ckpt
        --hic <hic_directory/>
  note that Hi-C data must be processed in npz format (can use [cool2npy.py](https://github.com/tanjimin/C.Origami/blob/main/src/corigami/preprocessing/cool2npy.py)), and stored in chromsome specific .npz files, analagous to the bigWigs described above.  
  
  To predict:  
      python3 predict.py 
      --<model_path>
      --mA <fire/>
      --mC <cpg/>
      --ctcf_fp <ctcf_footprint>
      --ctcf_dir <ctcf_motif_direction>



  
### Acknowledgements:  
The model architecture is inspired by, and largely a fork of, [C.Origami](https://github.com/tanjimin/C.Origami)  (Tan et. al. *Nature Biotechnology* 2023, https://www.nature.com/articles/s41587-022-01612-8).


### Cite:  
**Danilo Dubocanin, Anna Kalygina, J. Matthew Franklin, Cy Chittenden, Mitchell R Vollger, Shane Neph, Andrew B Stergachis, Nicolas Altemose. Integrating Single-Molecule Sequencing and Deep Learning to Predict Haplotype-Specific 3D Chromatin Organization in a Mendelian Condition. *bioRxiv*. 2025**




