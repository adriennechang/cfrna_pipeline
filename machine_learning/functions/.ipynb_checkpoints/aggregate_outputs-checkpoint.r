##------------------------------------
## SUMMARIZE TRAIN OUTPUTS
##------------------------------------

# Feature filtering used with every feature selection algorithm

AGG_TRAIN_OUTPUTS <- function(meta_data_train,
                              group1_name,
                              group2_name,
                              FILE_NAME,
                              PREFIX
                             ){
    
    #####
    ## input: 
    # - counts = count matrix, normalized cpm (Columns are genes, rows are samples)
    # - CORR_THRESH = max corerlation to not remove a gene
    # - LOW_THRESH = minimum sum of all counts for that gene to be removed
    ## output:
    # - counts.cor.var.cnts = filtered count matrix
    #####


    suppressMessages(library(pROC))
    suppressMessages(library(vioplot))
    suppressMessages(library(PRROC))
    
    FILE_PATH = paste0("./output/",FILE_NAME,"/",PREFIX,"/train/",PREFIX,"_",FILE_NAME) 

    ##------------------------------------
    # Plot ROC
    png(file=paste0(FILE_PATH,"_train_ROC.png"))

        rocobj1 <- plot.roc(meta_data_train$group, meta_data_train$classifier_score,
                            main = paste0("Training Set: ", group1_name," vs ",group2_name),
                            percent=TRUE,
                            ci=TRUE, boot.n=10000, boot.stratified=TRUE,
                            #	           show.thres=FALSE,
                            print.auc=TRUE,                  # compute AUC (of AUC by default)
                            #                   print.auc.pattern="%.2f",
                            print.thres="best",
                            print.thres.best.method="youden")

    dev.off()
    
    roc_auc = roc(meta_data_train$group ~ meta_data_train$classifier_score, plot = FALSE, print.auc = TRUE)
    data.frame("metric" = "auc", "value" = as.numeric(roc_auc$auc))   
    write.table( data.frame("metric" = "auc", "value" = as.numeric(roc_auc$auc)) ,
                paste0(FILE_PATH,"_train_roc-auc.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    youden_threshold <- coords(rocobj1, x="best", input="threshold", best.method="youden")$threshold

    ##------------------------------------
    # Plot Violin Plot
    png(file=paste0(FILE_PATH,"_train_violin.png"))

        vioplot(classifier_score~group, 
                data=meta_data_train, 
                outpch=NA, ylab="Classifier score", ylim=c(0,1), 
                main = paste0("Training Set: ", group1_name," vs ",group2_name),
                col=c('gray','dimgray'))

        abline(h=youden_threshold,col = "gray", lty = 4 )

        stripchart(classifier_score~group,
                   data=meta_data_train,  
                   vertical=TRUE, 
                   method="jitter", 
                   add=TRUE,
                   pch = 16, 
                   col=c(4,2))

    dev.off()

    ##------------------------------------
    # Plot PR Curve
    png(file=paste0(FILE_PATH,"_train_PR.png"))

    prcurve = pr.curve(scores.class0 = meta_data_train$classifier_score,
             weights.class0 = ifelse(meta_data_train$group == group1_name, 1,0), 
             curve = TRUE
            )

    plot(prcurve,color="black")

    dev.off()

    ##------------------------------------
    # Calculate confusion matrix and save output

    meta_data_train$class_prediction <-
      ifelse(meta_data_train$classifier_score > youden_threshold,
             group1_name,
             group2_name
      )

    results <- confusionMatrix(factor(meta_data_train$class_prediction), factor(meta_data_train$group))

    write.table(data.frame(as.table(results)),
                paste0(FILE_PATH,"_train_confusion-mtrx.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    write.table(data.frame(value = as.matrix(results,what="overall"))%>% rownames_to_column("stat"),
                paste0(FILE_PATH,"_train_confusion-overall.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    write.table(data.frame(value=as.matrix(results,what="classes")) %>% rownames_to_column("stat"),
                paste0(FILE_PATH,"_train_confusion-classes.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)
    
    
    return(youden_threshold)
    
    }



##------------------------------------
## SUMMARIZE TEST OUTPUTS
##------------------------------------

# Feature filtering used with every feature selection algorithm

AGG_TEST_OUTPUTS <- function(meta_data_test,
                             youden_threshold,
                             group1_name,
                             group2_name,
                             FILE_NAME,
                             PREFIX
                             ){
    
    #####
    ## input: 
    # - counts = count matrix, normalized cpm (Columns are genes, rows are samples)
    # - CORR_THRESH = max corerlation to not remove a gene
    # - LOW_THRESH = minimum sum of all counts for that gene to be removed
    ## output:
    # - counts.cor.var.cnts = filtered count matrix
    #####


    suppressMessages(library(pROC))
    suppressMessages(library(vioplot))
    suppressMessages(library(PRROC))
    
    FILE_PATH = paste0("./output/",FILE_NAME,"/",PREFIX,"/test/",PREFIX,"_",FILE_NAME)

   ##------------------------------------
    # Plot ROC
    png(file=paste0(FILE_PATH,"_test_ROC.png"))

        rocobj1 <- plot.roc(meta_data_test$group, meta_data_test$classifier_score,
                            main = paste0("Test Set: ", group1_name," vs ",group2_name),
                            percent=TRUE,
                            ci=TRUE, boot.n=10000, boot.stratified=TRUE,
                            #	           show.thres=FALSE,
                            print.auc=TRUE)

    dev.off()
    
    roc_auc = roc(meta_data_test$group ~ meta_data_test$classifier_score, plot = FALSE, print.auc = TRUE)
    data.frame("metric" = "auc", "value" = as.numeric(roc_auc$auc))   
    write.table( data.frame("metric" = "auc", "value" = as.numeric(roc_auc$auc)) ,
                paste0(FILE_PATH,"_test_roc-auc.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    ##------------------------------------
    # Plot Violin Plot
    png(file=paste0(FILE_PATH,"_test_violin.png"))

        vioplot(classifier_score~group, 
                data=meta_data_test, 
                outpch=NA, ylab="Classifier score", ylim=c(0,1), 
                main = paste0("Test Set: ", group1_name," vs ",group2_name),
                col=c('gray','dimgray'))

        abline(h=youden_threshold,col = "gray", lty = 4 )

        stripchart(classifier_score~group,
                   data=meta_data_test,  
                   vertical=TRUE, 
                   method="jitter", 
                   add=TRUE,
                   pch = 16, 
                   col=c(4,2))

    dev.off()

    ##------------------------------------
    # Plot PR Curve
#     png(file=paste0(FILE_PATH,"_test_PR.png"))

#     prcurve = pr.curve(scores.class0 = meta_data_test$classifier_score,
#              weights.class0 = ifelse(meta_data_test$group != group1_name, 0,1), 
#              curve = TRUE
#             )

#     plot(prcurve,color="black")

#     dev.off()

    ##------------------------------------
    # Calculate confusion matrix and save output

    meta_data_test$class_prediction <-
      ifelse(meta_data_test$classifier_score > youden_threshold,
             group1_name,
             group2_name
      )

    results <- confusionMatrix(factor(meta_data_test$class_prediction), factor(meta_data_test$group))

    write.table(data.frame(as.table(results)),
                paste0(FILE_PATH,"_test_confusion-mtrx.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    write.table(data.frame(value = as.matrix(results,what="overall"))%>% rownames_to_column("stat"),
                paste0(FILE_PATH,"_test_confusion-overall.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)

    write.table(data.frame(value=as.matrix(results,what="classes")) %>% rownames_to_column("stat"),
                paste0(FILE_PATH,"_test_confusion-classes.tsv"), 
                sep = "\t", 
                quote=F,
                row.names=F)
    
    
    }


##------------------------------------
## COMPARE OUTPUTS OF DATA
##------------------------------------

get_data <- function(MODEL,FILE_NAME,METHOD){
    
    path = paste0("./output/",FILE_NAME,"/",MODEL,"/",METHOD,"/",MODEL,"_",FILE_NAME,"_",METHOD)
    overall <- read.delim(paste0(path,"_confusion-overall.tsv"))
    classes <- read.delim(paste0(path,"_confusion-classes.tsv"))

    all <- rbind(overall,classes) 
    
    colnames(all)[2] <- paste0(MODEL,"_",METHOD)

    return(all)
}



COMPARE_OUTPUTS <- function(FILE_NAME){
    
    ## Get list of models
    all_models <- list.files(paste0("output/",FILE_NAME))
    all_models <- all_models[!grepl("\\.png|\\.tsv|\\.rds",all_models)]

    ## Get all training and test data
    all_train <- lapply(all_models,get_data,FILE_NAME,"train")
    df_train <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "stat", all.x = TRUE),
            all_train)

    all_test <- lapply(all_models,get_data,FILE_NAME,"test")
    df_test <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "stat", all.x = TRUE),
            all_test)

    df_all <- merge(df_train,df_test, by = "stat")


    df_train_m <- reshape2::melt(df_all,id.var = "stat") %>%
        mutate(method = gsub(".*_","",variable),
               model = gsub("_.*","",variable)) %>%
        select(!variable)

    df_train_m <- df_train_m %>% filter(stat %in% c("Accuracy",
                                                    "Kappa",
                                                    "Sensitivity",
                                                    "Specificity"))

    tmp_plt <- df_train_m %>%
        mutate(model = factor(model, levels = df_train_m %>% filter(method == "train" & stat == "Accuracy") %>% arrange(desc(value)) %>% pull(model))) %>%
        mutate(stat = factor(stat, levels = c("Specificity","Sensitivity","Kappa","Accuracy"))) %>%
        ggplot(aes(x=value, y=stat, color = method, shape = method)) + 
        geom_point(size = 1.5) + 
        facet_wrap(~model,ncol= 1, dir = "h")+
        scale_x_continuous(limits = c(0,1) , breaks = c(0,0.2,0.4,0.6,0.8,1))+
        theme_bw() +
        labs(title = paste0(group1_name," vs ", group2_name)) +
        theme(axis.title.y = element_blank(), 
              plot.title = element_text(hjust = 0.5))+
#         scale_size_manual(values=c(1,1.5))+
        scale_shape_manual(values = c(17,19))
    
    ggsave(paste0("output/",FILE_NAME,"/",FILE_NAME,"_dotplot.comparison.png"), tmp_plt)
    
    return(tmp_plt)
}


##################################################################
## Compare outputs accross analytes


get_aucs <- function(comp, model_df){
    
    res_df <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(res_df) <- c("algorithm","comparison","analyte","train","test")
    
    ## create dataframe
    if (model_df %>% filter(comparison == comp) %>% pull(cfrna) == 1){

        path = paste0("./output/cfrna_",comp,"/")

        all_algs <- list.files(path)
        all_algs <- all_algs[!grepl(".png|.rds",all_algs)]

        for (alg in all_algs){
            train_auc <- read.delim(paste0("./output/cfrna_",comp,"/",alg,"/train/",alg,"_cfrna_",comp,"_train_roc-auc.tsv"))[1,2]
            test_auc <- read.delim(paste0("./output/cfrna_",comp,"/",alg,"/test/",alg,"_cfrna_",comp,"_test_roc-auc.tsv"))[1,2]
            
            res_df[nrow(res_df)+1,] <- c(alg,comp,"cfrna",train_auc,test_auc)
        }
    }

    if (model_df %>% filter(comparison == comp) %>% pull(wbrna) == 1){

        path = paste0("./output/wbrna_",comp,"/")

        all_algs <- list.files(path)
        all_algs <- all_algs[!grepl(".png|.rds",all_algs)]

        for (alg in all_algs){
            train_auc <- read.delim(paste0("./output/wbrna_",comp,"/",alg,"/train/",alg,"_wbrna_",comp,"_train_roc-auc.tsv"))[1,2]
            test_auc <- read.delim(paste0("./output/wbrna_",comp,"/",alg,"/test/",alg,"_wbrna_",comp,"_test_roc-auc.tsv"))[1,2]

            res_df[nrow(res_df)+1,] <- c(alg,comp,"wbrna",train_auc,test_auc)
        }
    }

    return(res_df)
}



COMPARE_ANALYTES <- function(comp,model_df){
    
    res_df <-  get_aucs(comp,model_df) %>% reshape2::melt(id.vars = c("algorithm","comparison","analyte")) %>% rename(auc = value) %>% mutate(auc = as.numeric(auc)) %>% filter(algorithm != "LOG")

    res_df %>%
        mutate(algorithm = factor(algorithm, levels = res_df %>% 
                                                      group_by(algorithm) %>%
                                                      summarize(auc = mean(auc)) %>%
                                                      arrange(auc) %>%
                                                      pull(algorithm))) %>%
        mutate(var_analyte = paste0(analyte," ",variable)) %>%
        ggplot(aes(x = auc, y = algorithm, group = analyte, color = variable, shape = var_analyte)) +
        geom_point(size = 3,
                   stroke = 1.25,
                   position=position_jitterdodge(jitter.width = 0, jitter.height = 0, seed = 42)) +

        theme_bw() +
        theme(
            axis.title.x = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_blank(), 
            axis.text.y = element_text(size = 12, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 16)) +
        labs(title = comp)+

        scale_color_manual(breaks = c("train","test"), 
                           values = c("#00A1D5FF","#B24745FF"))+
        scale_shape_manual(breaks = c("cfrna test","cfrna train","wbrna test","wbrna train"), 
                            values = c(1,16,2,17)) +  

        geom_hline(yintercept = seq(1.5, length(unique(res_df$algorithm)) - 0.5), 
                   size = 0.5) +
        geom_hline(yintercept = seq(1, length(unique(res_df$algorithm))), 
                   size = 0.25, 
                   linetype = "dashed")


    ggsave(paste0("./output/",comp,"_analyte-comp.png"))

}