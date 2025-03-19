## FiberFold enables the prediction of Hi-C maps (2MB window, 10kb resolution) from a single assay, Fiber-seq (PacBio).  

### Dependencies

We've provided a .yml file to install the necessary environment, but if you are having issues you can install and use these dependencies.

    python 3.11.8
    pyfaidx
    pybigwig
    CUDA v12.1
    cuDNN v8.9
    torch==2.2.2+cu121
    torchmetrics==1.3.2
    torchvision==0.17.2
    pytorch-lightning==1.9.5

    
### First, generate bigWigs of FIRE, mCpG, CTCF footprinting, and CTCF motif direction.  
  [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) for generating CpG (default on PacBio output)  
  [fibertools](https://github.com/fiberseq/fibertools-rs) FIRE for generating accessible regulatory element track  
  [fibertools](https://github.com/fiberseq/fibertools-rs) footprint for generating CTCF footprinting information, CTCF Motif to use: MA0139.1.  
  [filterOverlapMotif.py](https://github.com/altemoselab/fiberFold/blob/main/misc/filterOverlapMotifs.py) will filter out overlapping motif based on score.  
  [encodeFP.py](https://github.com/altemoselab/fiberFold/blob/main/misc/encodeFP.py) will create a bigWig track from fibertools footprint output.  
  [encodeMotifDirection.py](https://github.com/altemoselab/fiberFold/blob/main/misc/encodeMotifDirection.py) will create a bigWig indicating CTCF motif direction.  

  
  All tracks should be divided by depth to be in range 0-1.
  To generate phased Hi-C maps, phase Fiber-seq reads prior to the steps above.
   
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
        
  To download model from our preprint:  

      wget https://us-west-1.console.aws.amazon.com/s3/object/altemoselab-share?region=us-west-1&bucketType=general&prefix=fiberFold/model.ckpt
      
  To predict:   
      
      python3 predict.py 
        --model <model_path> 
        --mA <fire/>
        --mC <cpg/>
        --ctcf_fp <ctcf_footprint>
        --ctcf_dir <ctcf_motif_direction>
        --ds <intervals_to_predict.bed>

To view predictions:
        
        with open('predict_out.pkl','rb') as handle :
            predicted_maps = pickle.load(handle)
        for x in predicted_maps:
            plt.imshow(x[0][0],cmap = 'magma_r', vmin = 0, vmax = 5)
            plt.show() 

![image](https://github.com/user-attachments/assets/7db23e82-644a-473c-bd6e-60b4177c0696)


  
### Acknowledgements:  
The model architecture is inspired by, and largely a fork of, [C.Origami](https://github.com/tanjimin/C.Origami)  (Tan et. al. *Nature Biotechnology* 2023, https://www.nature.com/articles/s41587-022-01612-8).


### Cite:  
**Danilo Dubocanin, Anna Kalygina, J. Matthew Franklin, Cy Chittenden, Mitchell R Vollger, Shane Neph, Andrew B Stergachis, Nicolas Altemose. Integrating Single-Molecule Sequencing and Deep Learning to Predict Haplotype-Specific 3D Chromatin Organization in a Mendelian Condition. *bioRxiv*. 2025**




