# ChromaFold

ChromaFold is a deep learning model that enables prediction of 3D contact maps from scATAC-seq data alone, by using pseudobulk chromatin accessibility and co-accessibility from scATAC-seq as well as predicted CTCF motif tracks as input features. 

<img width="927" alt="model" src="https://github.com/viannegao/ChromaFold/assets/111778845/7eb5bdad-7547-4bc8-aab2-1db4df47ae1a">

<br/>
<br/>

## Requirements
 
- General
  - python=3.8
  - pytorch=1.11
  - numpy=1.21
  - pandas=1.4
  - scipy=1.7


- Visualization
  - coolbox=0.3
  - matplotlib=3.2
  - seaborn=0.11
  - tabix=1.11

You can create a conda environment using the provided .yml file:
```
conda env create -f chromafold.yml
```
 <br/>

## Data Preprocessing

**Raw data preparation**

Sample raw and processed input data can be downloaded from https://drive.google.com/drive/folders/1p6dulb2z51NF_WA6RnAG4hHuUaKfFPrR?usp=sharing

a) Input data preparation
  - Prepare CTCF motif data: CTCF motif data are extracted from the [CTCF introduction](https://bioconductor.org/packages/release/data/annotation/vignettes/CTCF/inst/doc/CTCF.html) from R package [AnnotationHub](https://www.bioconductor.org/packages/release/bioc/html/AnnotationHub.html). R scripts for generating motif of hg38 and mm10 can be found at [process_input/ctcf_motif](https://github.com/viannegao/ChromaFold/tree/main/process%20input/ctcf_motif). We also provide ready-to-use CTCF motif score for hg38, hg19, mm10 in the [google drive](https://drive.google.com/drive/folders/1TlfXGix2U-K1UIrOYv8gWGIiSx10dP3M).
  - Prepare scATAC data for inference: please refer to the full instructions at [preprocessing_pipeline](https://github.com/viannegao/ChromaFold/tree/main/preprocessing_pipeline).

b) Target data preparation

  - Example raw Hi-C file for IMR-90 can be downloaded from ENCODE (https://www.encodeproject.org/files/ENCFF843MZF/@@download/ENCFF843MZF.hic).
  - Prepare normalized Hi-C library for target: please refer to the full instructions at [process input/hic_normalization](https://github.com/viannegao/ChromaFold/tree/main/process%20input/hic_normalization).

**Integration for training**

  - Prepare Hi-C data for training
    
    - Run "process input/Process Input - Hi-C.ipynb". The juicer tools jar file can be downloaded from https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar . If the juicer tool doesn't match your java system, please refer to [earlier versions](https://github.com/aidenlab/juicer/wiki/Download) of the juicer tools. 
    
   

 <br/>

## Inference
 <br/>


**1. Run inference on germinal center B cell with ChromaFold**
  - Run inference on full chromosome without offset
  ```
  python script/inference.py --data-path ./data/processed_input/ -ct gcb --model-path ./checkpoint/chromafold_CTCFmotif.pth.tar -chrom 16 -offset -2000000 --genome mm10
  ```
  - Run inference only on regions with complete input information
  ```
  python script/inference.py --data-path ./data/processed_input/ -ct gcb --model-path ./checkpoint/chromafold_CTCFmotif.pth.tar -chrom 16 -offset 0 --genome mm10
  ```
  
## Training
 <br/>

**1. Training on 3 cell types**
  - Train model without co-accessibility component
  ```
  python script/train_bulkOnly.py --data-path ./data/processed_input/ -ct gm12878_hg38 umb_endo imr90
  ```
  - Train model with co-accessibility
  ```
  python script/train.py --data-path ./data/processed_input/ -ct gm12878_hg38 umb_endo imr90 --mod-name bothInput
  ```
  - Train deterministic model for full reproducibility
  ```
  python script/train_bulkOnly.py --data-path ./data/processed_input/ -ct gm12878_hg38 umb_endo imr90 --deterministic --mod-name deterministic
  ```
  
 <br/>
 
**2. Training on 1 cell type**
  - Train model on HUVEC without co-accessibility component
  ```
  python script/train_bulkOnly.py --data-path ./data/processed_input/ -ct umb_endo
  ```

