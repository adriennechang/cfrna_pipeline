##------------------------------------
## LOGISTIC REGRESSION
##------------------------------------

GLM_TRAIN <- function(counts,meta_data){
    
    set.seed(42)
    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    glm_fit <- train(y~., 
                    data=counts, 
                    method='glm',
                    preProcess = c("center","scale"),
                     weights = model_weights,
                    family = "binomial")
    
    return(glm_fit)

}



##------------------------------------ X
## RANDOM FOREST
##------------------------------------

RF_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)

    metric <- "Accuracy"
    mtry <- sqrt(ncol(counts)-1)
    tunegrid <- expand.grid(.mtry = (5:50))
    
    counts$y <- meta_data$group

    rf_fit <- train(y~., 
                          data=counts, 
                          method='rf', 
                          metric='Accuracy', 
                          tuneGrid = tunegrid, 
                          trControl=control)
    
    
    return(rf_fit)

}


##------------------------------------ X
## SUPPORT VECTOR MACHINE LINEAR
##------------------------------------

SVMLIN_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5,classProbs = TRUE)
    
    tunegrid <- expand.grid(C = seq(0.0001, 2, length = 20))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    svm_fit <- train(y~., 
                    data=counts, 
                    method='svmLinear', 
                    metric = metric,
                    trControl=control,
                     tuneGrid = tunegrid,
                    preProcess = c("center","scale"))
    
    return(svm_fit)

}


##------------------------------------ X
## Naive Bayes
##------------------------------------

NB_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(fL = 0:5,
                            usekernel = c(TRUE,FALSE),
                            adjust = seq(0, 5, by = 1)
                            )

    counts$y <- factor(meta_data$group)

    nb_fit <- train(y~., 
                    data=counts, 
                    method='nb', 
                    metric = metric,
                    preProcess = c("center","scale"),
                    trControl=control,
                   tuneGrid = tunegrid)
    
    return(nb_fit)

}


##------------------------------------
## LDA
##------------------------------------

LDA_TRAIN <- function(counts,meta_data){
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    counts$y <- factor(meta_data$group)

    lda_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    preProcess = c("center","scale"),
                    method='lda', 
                    trControl=control)
    
    return(lda_fit)

}


##------------------------------------
## KNN
##------------------------------------

KNN_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)
    
    counts$y <- factor(meta_data$group)

    knn_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='knn', 
                    trControl=control,
                    preProcess = c("center","scale"),
                    tuneLength = 20                      # test different K's
                    )
    
    return(knn_fit)

}

##------------------------------------
## GLMNETRidge
##------------------------------------

GLMNETRidge_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    tunegrid <- expand.grid(alpha = 0,
                           lambda = seq(0.0001, 1, length = 100))
    
    counts$y <- factor(meta_data$group)

    glmnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='glmnet',
                    preProcess = c("center","scale"),
                    trControl=control)
    
    return(glmnet_fit)
    
}


##------------------------------------
## GLMNETLasso
##------------------------------------

GLMNETLasso_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    tunegrid <- expand.grid(alpha = 1,
                           lambda = seq(0.0001, 1, length = 100))
    
    counts$y <- factor(meta_data$group)

    glmnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='glmnet',
                    preProcess = c("center","scale"),
                    trControl=control)
    
    return(glmnet_fit)
    
}


##------------------------------------
## RPART
##------------------------------------

RPART_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    rpart_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='rpart', 
#                     preProcess = c("center","scale"),
                    trControl=control,
                    weights = model_weights,
                    tuneLength = 20)
    
    return(rpart_fit)

}


##------------------------------------
## XGBTREE
##------------------------------------

XGBTREE_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)

    
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    xgbtree_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='xgbTree', 
                    preProcess = c("center","scale"),
                    weights = model_weights,
                    trControl=control)
    
    return(xgbtree_fit)

}

##------------------------------------
## NNET
##------------------------------------

NNET_TRAIN <- function(counts,meta_data){
    
    
    set.seed(42)
    
    metric <- "Accuracy"

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                            size = c(3, 5, 10, 20))

    
    counts$y <- factor(meta_data$group)
    
   model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    nnet_fit <- train(y~., 
                    data=counts, 
                    metric = metric,
                    method='nnet', 
                    preProcess = c("center","scale"),
                    trControl=control,
                    tuneGrid = tunegrid,
                      weights = model_weights,
                     trace = FALSE)
    
    return(nnet_fit)

}

##------------------------------------
## SVM RADIAL
##------------------------------------

SVMRAD_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5,classProbs = TRUE)
    
    tunegrid <- expand.grid(C = seq(0, 2, length = 20),
                           sigma = c(.01, .015, 0.2))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    svmrad_fit <- train(y~., 
                    data=counts, 
                    method='svmRadial', 
                    metric = metric,
                    trControl=control,
		    tuneGrid=tunegrid, ## AC ADD
                    preProcess = c("center","scale"))
    
    return(svmrad_fit)

}

##------------------------------------
## PAM
##------------------------------------

PAM_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    pam_fit <- train(y~., 
                    data=counts, 
                    method='pam', 
                    metric = metric,
                    trControl=control,
                    tuneLength = 20,
                    preProcess = c("center","scale"))
    
    return(pam_fit)

}

##------------------------------------
## ADABOOST
##------------------------------------

ADABOOST_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(nlter = seq(0, 2, length = 20),
                           method = c(.01, .015, 0.2))
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)

    aba_fit <- train(y~., 
                    data = counts, 
                    method = 'adaboost', 
                    metric = metric,
                    trControl = control,
                    tuneLength = 30,
                    preProcess = c("center","scale"))
    
    return(aba_fit)

}

##------------------------------------
## C5.0
##------------------------------------

C5_TRAIN <- function(counts,meta_data){
    
    set.seed(42)

    control <- trainControl(method="cv",number=5)
    
    tunegrid <- expand.grid(winnow = c(TRUE,FALSE), 
                            trials=c(1,5,10,15,20),
                            model="tree" )
    
    
    metric <- "Accuracy"
    counts$y <- factor(meta_data$group)
    
    model_weights <- ifelse(counts$y == names(table(counts$y)[1]),
                    (1/table(counts$y)[1]) * 0.5,
                    (1/table(counts$y)[2]) * 0.5)

    c5_fit <- train(y~., 
                    data = counts, 
                    method = 'C5.0', 
                    metric = metric,
                    trControl = control,
                    tuneGrid = tunegrid,
                    weigths = model_weights,
                    preProcess = c("center","scale"))
    
    return(c5_fit)

}


##------------------------------------ X
## Extra Trees
##------------------------------------

EXTRATREES_TRAIN <- function(counts,meta_data){
    
    set.seed(42)
  require(ranger)
    control <- trainControl(method="cv",number=5, 
                            classProbs = TRUE)
        
    tunegrid <- expand.grid(mtry = seq(5, ncol(counts), by = 5),
                           min.node.size = seq(2,10),
                            splitrule = "extratrees"
                           )

    metric <- "Accuracy"
    
    counts$y <- meta_data$group

    rf_fit <- train(y~., 
                          data=counts, 
                          method='ranger', 
                          metric='Accuracy', 
                          tuneGrid = tunegrid, 
                          trControl=control,
#                           preProc = c("center", "scale"),
                        na.action = na.pass
                          )
    
    return(rf_fit)
}

