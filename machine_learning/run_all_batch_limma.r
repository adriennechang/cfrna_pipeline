suppressMessages(library(tidyverse))
suppressMessages(library(caret))
source("./functions/main_training_tof.r")
source("./functions/training_models.r")
source("./functions/aggregate_outputs.r")
source("./functions/feature_selection_methods_batch_limma.r")

##------------------------------------
## Set paths to necessary files
##------------------------------------


count_matrix_path = "../files_for_manuscript/count_matrix.txt"
meta_data_path = "../files_for_manuscript/meta_data.txt"
id_column <- "Cornell_ID"                                             ## Column which matches to rownames in count matrix
group_column = "Microbiologic.reference.standard"
## References for gene filtering
gene_blacklist_path <- "../files_for_manuscript/genes_to_exlude.tsv"
gene_biotype_key_path <- "../files_for_manuscript/gencode.biotype.name.key.tsv"

##------------------------------------
## Set parameters
##------------------------------------

NGENES = 150                                                      ## Number of top differentially abundant genes to use
TEST_TRAIN_SPLIT = 0.7                                                 ## how to split data for train and test


##------------------------------------
## Run comparisons
##------------------------------------



FILE_NAME = "test_paper"

group1_name = "Positive"
group2_name = "Negative"



stats <- MAIN(group1_name,
     group2_name,
     FILE_NAME,
     count_matrix_path,
     meta_data_path,
     id_column,
     group_column, 
     gene_blacklist_path,
     gene_biotype_key_path,
     NGENES,
     TEST_TRAIN_SPLIT)

comp_stats <- COMPARE_OUTPUTS(FILE_NAME)

