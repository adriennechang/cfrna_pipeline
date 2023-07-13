# Circulating Cell-Free RNA in Blood as a Host Response Biomarker for the Detection of Tuberculosis

This repository contains the pipeline and analysis scripts used to evaluate TB markers in plasma cfRNA.

## Overview
To run the cfRNA pipeline the following must be done:
1. Create the conda environment
2. Copy the pipeline
3. Execute test run

### Create the conda environment
1. Create the environment from the environment yml file.
```
conda env create -f cfrna_paper_env.yml
```
2. Activate the new environment.
```
conda activate cfrna_pipeline
```


### Copy the pipeline and run a test:
Estimated time: 20 minutes 
1. Clone the repository
```
git clone https://github.com/adriennechang/cfrna_pipeline
```
2. Change into the pipeline directory and make a sequencing prep table.
3. Run snakemake
```
snakemake --cores 2 -k
```

## Machine learning pipeline
The machine learning pipeline can be performing the following steps:
```
cd machine_learning/
Rscript run_all_batch_limma.r
```

## Analysis
The expected outputs from the main pipeline and machine learning pipeline can be found in "files_for_manuscript". The following software and versions were used to generate the analysis scripts, which are provided for each figure:  
  
R version 4.1.3 (2022-03-10)  
Platform: x86_64-conda-linux-gnu (64-bit)  
Running under: Rocky Linux 9.0 (Blue Onyx)  

- reshape2_1.4.4
- ggplot2_3.3.6 
- qualpalr_0.4.3
- dplyr_1.0.10 
- tidyverse_1.3.2
- DESeq2_1.34.0 
- data.table_1.14.4
- pheatmap_1.0.12
- caret_6.0-93
- pROC_1.18.0
- edgeR_3.36.0
- limma_3.50.3 
- magrittr_2.0.3
- stringr_1.4.1 
- ggpubr_0.4.0

