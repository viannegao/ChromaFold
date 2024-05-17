## Data Preprocessing for ChromaFold 

All processing steps have been summarized in [`fragment_to_input.sh`](https://github.com/viannegao/ChromaFold/blob/main/preprocessing_pipeline/fragment_to_input.sh) file. Below is a more detailed explanation of the steps. 

### **Step 0**. Filter and merge fragement files

- **Input:** single-cell ATAC-seq fragment file and barcode file. Example data: [google drive](https://drive.google.com/drive/folders/1ZDwumdoC-9lqsVEHUeBs4euN8xcBelhu?usp=sharing)
- **Script:** [`fragment_celltype_merge.py`](https://github.com/viannegao/ChromaFold/blob/main/preprocessing_pipeline/fragment_celltype_merge.py)

Filter and merge fragments of a specific cell type from multiple fragment datasets. In the example data from google drive, since we only have a single scATAC-seq fragment dataset, we will skip this step. 

### **Step 1**. Run [ArchR](https://www.archrproject.com/bookdown/co-accessibility-with-archr.html) to calculate co-accessibility

- **Input:** scATAC-seq fragment file
- **Package:** [ArchR](https://www.archrproject.com/bookdown/co-accessibility-with-archr.html)



