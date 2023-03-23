# Machine Learning Testing
# Conor Loy 2022-08-18 
# This code will perform feature selection and training/cross-validation


MAIN <- function(group1_name,
                 group2_name,
                 FILE_NAME,
                 count_matrix_path,
                 meta_data_path,
                 id_column,
                 group_column,
                 gene_blacklist_path,
                 gene_biotype_key_path,
                 NGENES,
                 TEST_TRAIN_SPLIT
                ){
    
    cat("-------------------------------------\n")
    cat(FILE_NAME)
    cat("\n\n")
    
    
    if (!dir.exists(paste0("output/",FILE_NAME))){
    system(paste0("mkdir output/",FILE_NAME))}                           ## Make output folder
    
    ##------------------------------------
    ## READ IN META DATA AND COUNT MATRIX
    ##------------------------------------
    cat("--> READING IN DATA AND SPLITTING\n")
    read_data_list <- READ_DATA(group1_name,
                 group2_name,
                 count_matrix_path,
                 meta_data_path,
                 id_column,
                 group_column,
                 gene_blacklist_path)
    
    count_matrix <- read_data_list[['count_matrix']]
    meta_data <- read_data_list[['meta_data']]
    

    if (!all(colnames(count_matrix) == meta_data$Cornell_ID)){
        print("ISSUE")
        count_matrix = data.frame()
    }
    
    ##------------------------------------
    ## SPLIT DATA INTO TRAIN AND TEST
    ##------------------------------------
    split_data_list <- SPLIT_DATA(count_matrix,
                                  meta_data,
                                  group1_name,
                                 TEST_TRAIN_SPLIT)
    
    meta_data_train <- split_data_list[['meta_data_train']]
    meta_data_test <- split_data_list[['meta_data_test']]
    count_matrix_train <- split_data_list[['count_matrix_train']]
    count_matrix_test <- split_data_list[['count_matrix_test']]
    
    ##------------------------------------
    ## PERFORM FEATURE SELECTION
    ##------------------------------------
    cat("--> PERFORMING FEATURE SELECTION\n")
    feature_selection_list <- FEATURE_SELECTION_NORMALIZATION(count_matrix_train,
                                                count_matrix_test,
                                                meta_data_train,
                                                group1_name,
                                                group2_name,
                                                NGENES)
    
    count_matrix_train_sig = feature_selection_list[["count_matrix_train_sig"]]
    count_matrix_test_sig = feature_selection_list[["count_matrix_test_sig"]]
    
    ##------------------------------------
    ## TRAIN ALL MODELS
    ##------------------------------------
    cat("--> TRAINING MODELS\n")
    suppressWarnings(suppressMessages(TRAIN_ALL_MODELS(count_matrix_train_sig,
                     count_matrix_test_sig,
                     meta_data_train,
                     meta_data_test,
                     group1_name,
                     group2_name,
                     FILE_NAME)))
    
    
    saveRDS(list(
            "meta_data_train" = meta_data_train,
            "meta_data_test" = meta_data_test,
            "count_matrix_train" = count_matrix_train,
            "count_matrix_test" = count_matrix_test,
            "count_matrix_train_sig" = count_matrix_train_sig,
            "count_matrix_test_sig" = count_matrix_test_sig
            ),
                paste0("output/",FILE_NAME,"/",FILE_NAME,"_data.rds")
           )
    
    cat("-------------------------------------\n")
}


##------------------------------------------------------------------------------------------------------------------------------------------------------------
READ_DATA <- function(group1_name,
                 group2_name,
                 count_matrix_path,
                 meta_data_path,
                 id_column,
                 group_column,
                 gene_blacklist_path){
    
    ##------------------------------------
    ## Read in meta data and count matrix, remove blacklisted genes, and synchronize order
    ##------------------------------------
    suppressMessages(library(tidyverse))

    ##------------------------------------
    ## Meta Data
	meta_data = read.delim(meta_data_path)
   meta_data[[group_column]] = gsub(meta_data[[group_column]], pattern="TB ", replacement="")

    meta_data$group <- factor(meta_data[[group_column]], levels = c(group1_name,group2_name))
	meta_data = meta_data[!is.na(meta_data$group),]
 #meta_data = meta_data[meta_data$HIV.status=="Negative",]
    table(meta_data$group)

    ##------------------------------------
    ## Count Matrix
    gene_blacklist <-  read.delim(gene_blacklist_path,col.names = c("type,","ENSMBL","gene_symbol"))   ## read in gene list to exclude
    cols = intersect(colnames(read.table(count_matrix_path, row.names=c(1))), meta_data$Cornell_ID)
    
    count_matrix <- read.table(count_matrix_path, row.names=c(1)) %>%
                    data.frame()

    gene.ids <- gsub("\\..*","",rownames(count_matrix))                                                ## clean up any weird gene names

    exclude.idx <- gene.ids %in% gene_blacklist[,2]                                                       

    count_matrix = count_matrix[!exclude.idx,]                                                         ## exclude genes

    meta_data = meta_data[meta_data$Cornell_ID %in% cols,]
    count_matrix <- count_matrix[,meta_data$Cornell_ID]                                                 ## Synchronize count and meta data frames
    
    return(list("meta_data" = meta_data, "count_matrix" = count_matrix))
    
}


##------------------------------------------------------------------------------------------------------------------------------------------------------------
SPLIT_DATA <- function(count_matrix,
                                  meta_data,
                                  group1_name,
                                 TEST_TRAIN_SPLIT){
    ##------------------------------------
    ## Partition Data into train and test sets
    ##------------------------------------
    #suppressMessages(library(caret))
    #suppressMessages(library(tibble))
    # set.seed(42)

    # trainIndex <- createDataPartition(meta_data$group, 
    #                                   p = TEST_TRAIN_SPLIT, 
    #                                   list = FALSE, 
    #                                   times = 1)

    # ## training data
    # meta_data_train <- meta_data[trainIndex,]
    # count_matrix_train <- count_matrix[,trainIndex]

    # ## test data
    # meta_data_test <- meta_data[-trainIndex,]
    # count_matrix_test <- count_matrix[,-trainIndex]


    set.seed(42)

    ## Slice samples by sample group
    meta_data_train <- meta_data %>% slice_sample(prop=0.7)

  
    ## create test data based on left overs
    meta_data_test <- meta_data %>% filter(!(Cornell_ID %in% meta_data_train$Cornell_ID))
    
    ## Subset count matrixs
    count_matrix_train <- count_matrix[,meta_data_train$Cornell_ID]
    count_matrix_test <- count_matrix[,meta_data_test$Cornell_ID]
    
    ## rename rows
    rownames(count_matrix_train) <- gsub("\\_.*","",rownames(count_matrix_train))
    rownames(count_matrix_test) <- gsub("\\_.*","",rownames(count_matrix_test))
    
    return(list('meta_data_train' = meta_data_train,
   'meta_data_test' = meta_data_test,
   'count_matrix_train' = count_matrix_train,
    'count_matrix_test' = count_matrix_test))
    }



##------------------------------------------------------------------------------------------------------------------------------------------------------------
meaure_performance <- function(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               TYPE){
        
    ###
    
    if (!dir.exists(paste0(FILE_NAME,"/",output_subfolder))){
        system(paste0("mkdir ",paste0("output/",FILE_NAME,"/",output_subfolder)))
        system(paste0("mkdir ",paste0("output/",FILE_NAME,"/",output_subfolder,"/train")))
        system(paste0("mkdir ",paste0("output/",FILE_NAME,"/",output_subfolder,"/test")))
    }    
    saveRDS(MODEL_FIT,paste0("output/",FILE_NAME,"/",output_subfolder,"/",output_subfolder,".model"))

    ##------------------------------------
    # Make predictions on training data
    tpred = predict(MODEL_FIT, newdata = count_matrix_train_sig, type=TYPE)

    meta_data_train$classifier_score <- as.numeric(tpred[,1])


    youden_threshold <- AGG_TRAIN_OUTPUTS(meta_data_train,
                      group1_name,
                      group2_name,
                      FILE_NAME,
                      output_subfolder)

    ##------------------------------------
    # TEST
    ##------------------------------------
    
    tpred = predict(MODEL_FIT, newdata = count_matrix_test_sig, type=TYPE)

    meta_data_test$classifier_score <- as.numeric(tpred[,1])

    AGG_TEST_OUTPUTS(meta_data_test,
                     youden_threshold,
                     group1_name,
                     group2_name,
                      FILE_NAME,
                      output_subfolder)    
}


##------------------------------------------------------------------------------------------------------------------------------------------------------------

TRAIN_ALL_MODELS <- function(count_matrix_train_sig,
                     count_matrix_test_sig,
                     meta_data_train,
                     meta_data_test,
                     group1_name,
                     group2_name,
                     FILE_NAME){
        
    all_models <- list()
    
##----------------------------------------------------------------------------
## RF
##----------------------------------------------------------------------------
    
    output_subfolder <- "RF"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.rf = RF_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.rf
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## SVM Linear
##----------------------------------------------------------------------------
    output_subfolder <- "SVMLin"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.svmlin = SVMLIN_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.svmlin
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## NB
##----------------------------------------------------------------------------
    output_subfolder <- "NB"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.nb = suppressWarnings(NB_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.nb
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## LDA
##----------------------------------------------------------------------------
    output_subfolder <- "LDA"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.lda = suppressWarnings(LDA_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.lda
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
    
##----------------------------------------------------------------------------
## KNN
##----------------------------------------------------------------------------
    output_subfolder <- "KNN"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.knn = KNN_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.knn
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
    
##----------------------------------------------------------------------------
## GLM]
##----------------------------------------------------------------------------
    output_subfolder <- "GLM"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.glm = GLM_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.glm
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

#     all_models[[output_subfolder]] <- MODEL_FIT
      
        
##----------------------------------------------------------------------------
## GLMNET RIDGE
##----------------------------------------------------------------------------
    output_subfolder <- "GLMNETRidge"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.glmnet = GLMNETRidge_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.glmnet
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## GLMNET LASSO
##----------------------------------------------------------------------------
    output_subfolder <- "GLMNETLasso"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.glmnet = GLMNETLasso_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.glmnet
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
    
##----------------------------------------------------------------------------
## RPART NET
##----------------------------------------------------------------------------
    output_subfolder <- "RPART"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.rpart = RPART_TRAIN(count_matrix_train_sig,meta_data_train)
    MODEL_FIT <- fit.rpart
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
    
##----------------------------------------------------------------------------
## EXTRATREES_TRAIN
##----------------------------------------------------------------------------
    output_subfolder <- "EXTRATREES"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.extratrees = suppressWarnings(EXTRATREES_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.extratrees
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT

##----------------------------------------------------------------------------
## NNET
##----------------------------------------------------------------------------
    output_subfolder <- "NNET"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.nnet = suppressWarnings(NNET_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.nnet
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
 
    
##----------------------------------------------------------------------------
## SVM Radial
##----------------------------------------------------------------------------
    output_subfolder <- "SVMRAD"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.svmrad = suppressWarnings(SVMRAD_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.svmrad
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## PAM
##----------------------------------------------------------------------------
    output_subfolder <- "PAM"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.pam = suppressWarnings(PAM_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.pam
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## ADABOOST
##----------------------------------------------------------------------------
#     output_subfolder <- "ADABOOST"
#     cat(paste0("    |--> ",output_subfolder,"\n"))
#     fit.pam = suppressWarnings(ADABOOST_TRAIN(count_matrix_train_sig,meta_data_train))
#     MODEL_FIT <- fit.pam
    
#     meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
#                                count_matrix_train_sig,count_matrix_test_sig,
#                                meta_data_train,meta_data_test,
#                                group1_name,group2_name,
#                                "prob")

#     all_models[[output_subfolder]] <- MODEL_FIT
    
    
##----------------------------------------------------------------------------
## C5
##----------------------------------------------------------------------------
    output_subfolder <- "C5"
    cat(paste0("    |--> ",output_subfolder,"\n"))
    fit.c5 = suppressWarnings(C5_TRAIN(count_matrix_train_sig,meta_data_train))
    MODEL_FIT <- fit.c5
    
    meaure_performance(MODEL_FIT,FILE_NAME,output_subfolder,
                               count_matrix_train_sig,count_matrix_test_sig,
                               meta_data_train,meta_data_test,
                               group1_name,group2_name,
                               "prob")

    all_models[[output_subfolder]] <- MODEL_FIT
    
##----------------------------------------------------------------------------
## Save Final results
##----------------------------------------------------------------------------
    
    ## Get accuracy and kappa statistics
    results <- resamples(all_models)
    result_sum <- summary(results)
    accuracy <- data.frame(result_sum$statistics$Accuracy)
    accuracy$model <- rownames(accuracy)
    accuracy$stat = "Accuracy"

    kappa <- data.frame(result_sum$statistics$Kappa)
    kappa$model <- rownames(kappa)
    kappa$stat = "Kappa"

    stats = rbind(accuracy,kappa)
    # calculate error 95%
    values <- results$values
    calc_error_95 <- function(X){ qnorm(0.95) * sd(X)/sqrt(length(X))}
    errors <- apply(values,2,calc_error_95) %>%
                data.frame() %>% rownames_to_column("model_stat") %>% rename(. = "error_95") %>%
                mutate(model = gsub("\\~.*","",model_stat), stat = gsub(".*\\~","",model_stat)) %>% 
                select(!model_stat) %>% filter(model!="Resample")

    # merge
    stats <- merge(stats, errors, by = c("model","stat"))
    
    write.table(stats,paste0("output/",FILE_NAME,"/",output_subfolder,"/",output_subfolder,".stats.tsv"),
                sep = "\t", row.names=F,quote=F)
    
}

