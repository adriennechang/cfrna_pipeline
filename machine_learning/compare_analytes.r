suppressMessages(library(tidyverse))
source("./functions/main_training.r")

#####################################################################################
#####################################################################################

source("./functions/aggregate_outputs.r")

## COMPARE ALL OUTPUTS
all_comparisons <- list.files("./output")
groups <- unique(gsub("cfrna_|wbrna_","",all_comparisons))

print(groups)

model_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(model_df) <- c("comparison","cfrna","wbrna")

for (g in groups){
    cfrna <- ifelse(file.exists(paste0("./output/cfrna_",g)), 1, 0)
    wbrna <- ifelse(file.exists(paste0("./output/wbrna_",g)), 1, 0)
    
    model_df[nrow(model_df)+1,] <- c(g,cfrna,wbrna)
}

head(model_df)



lapply(model_df$comparison, COMPARE_ANALYTES, model_df)