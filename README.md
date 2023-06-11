# ChromaFold

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
