## Hi-C normalization 

This page is for preparing [HiC-DC+](https://github.com/mervesa/HiCDCPlus) normalized Hi-C libraries for training. 

### **Step 0**. Generate HiC-DC+ normalized z-score files

- **Input:** HiC-DC+ requires assembly and resolution specific feature files for regression. Ready-to-use genomic features have been uploaded to [feature files](https://drive.google.com/drive/folders/1084P15MIrYeS13_ynpx2fbu11HTDszOt?usp=sharing). Please download them for normalization. 
- **Software and tools:** please download [Straw](https://github.com/aidenlab/straw) for data extraction from `.hic` format. A ready-to-use `straw.cpp` can also be downloaded from [google drive](https://drive.google.com/drive/folders/11BSDjUA4fb9uLnAqSFXjKZWQbzVadSIn?usp=sharing).
- **Script: ** [`hicdcplus_run.R`]()

```
Rscript /chromafold/process input/hic_normalization/hicdcplus_run.R \
10000 \
"/scripts/hicdc_workflow/hg38_MboI_10kb_features.rds" \
"/data/HiC/imr90/hicdc_out/" \
"Hsapiens" \
"hg38" \
"/data/HiC/imr90/ENCFF281ILS.hic" \
"scripts/hicdc_workflow/straw.cpp" \
"/scripts/hicdc_workflow/juicer_tools.1.9.9_jcuda.0.8.jar" 
```