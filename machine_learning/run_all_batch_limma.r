suppressMessages(library(tidyverse))
suppressMessages(library(caret))
source("./functions/main_training_tof.r")
source("./functions/training_models.r")
source("./functions/aggregate_outputs.r")
source("./functions/feature_selection_methods_batch_limma.r")

##------------------------------------
## Set paths to necessary files
##------------------------------------

## Count matrix

#count_matrix_path = "../takara_human_V2/resequence/output/GHL/GHL_feature_counts.tsv"
## Meta data
count_matrix_path = "../takara_human_V2/resequence/output/GHL_R2D2_feature_counts.tsv"
meta_data_path = "../analysis/GHL_R2D2_metadata_thresholds.txt"
#meta_data_path <- "../analysis/GHL_R2D2_metadata.txt"
#qc = '../takara_human_V2/resequence/output/GHL/GHL_mapping_stats.tsv'
qc = '../takara_human_V2/resequence/output/GHL_R2D2_P3_mapping_stats_QC.tsv'
id_column <- "Cornell_ID"                                             ## Column which matches to rownames in count matrix
#group_column <- "Semiquant.Xpert.result.2"  
group_column = "Microbiologic.reference.standard"
## References for gene filtering
gene_blacklist_path <- "/workdir/cfrna/references/human/hg38/genes_to_exlude.tsv"
gene_biotype_key_path <- "/workdir/cfrna/references/human/hg38/gencode.biotype.name.key.tsv"

##------------------------------------
## Set parameters
##------------------------------------

NGENES = 150                                                      ## Number of top differentially abundant genes to use
TEST_TRAIN_SPLIT = 0.7                                                 ## how to split data for train and test


##------------------------------------
## Run comparisons
##------------------------------------



FILE_NAME = "TB_P3_combined_batch_limma_tof"

group1_name = "Positive"
group2_name = "Negative"



stats <- MAIN(group1_name,
     group2_name,
     FILE_NAME,
     count_matrix_path,
     meta_data_path,qc,
     id_column,
     group_column, 
     gene_blacklist_path,
     gene_biotype_key_path,
     NGENES,
     TEST_TRAIN_SPLIT)

comp_stats <- COMPARE_OUTPUTS(FILE_NAME)

