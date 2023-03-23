# Circulating Cell-Free RNA in Blood as a Host Response Biomarker for the Detection of Tuberculosis

This repository contains the pipeline and analysis scripts used to evaluate TB markers in plasma cfRNA.

## Overview
To run the cfRNA pipeline the following must be done:
1. Install miniconda3
2. Create the conda environment
3. Copy the pipeline
4. Execute test run

### Install miniconda3:
Estimated time: 30 minutes

1. cd to your home directory
```
cd ~
```
2. Download and run the installation script
```
wget
https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
 ```
### Create the conda environment
Estimated time: a few hours
1. Activate conda
2. Create and activate a new conda environment
```
conda create -c bioconda -n cfrna_pipeline pip git R=4.1.0 python=3.7.10 fastqc=0.11.8 samtools=1.15.1 snakemake=7.7.0 umi_tools=1.1.2 subread=2.0.1 kraken2=2.1.0 bracken=2.7 r-data.table=1.14.2 ete3=3.1.2 r-devtools r-tidyverse=1.3.1
krakenuniq=0.7.3?
conda activate cfrna_pipeline
``` 
3. Install extra python packages using pip
```
pip install argparse==1.1
```
5. Install extra R packages in a R interactive session
```
R
library(“devtools”)
install_github("Danko-Lab/BayesPrism/BayesPrism")
q()
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
The following software and versions were used to generate the analysis scripts:  
  
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

