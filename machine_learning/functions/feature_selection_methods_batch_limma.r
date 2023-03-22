# Feature selection algorithms used in "main_training.r"
# Conor Loy 2022-08-18 

##------------------------------------------------------------------------------------------------------------------------------------------------------------
FEATURE_SELECTION_NORMALIZATION <- function(count_matrix_train,
                                                count_matrix_test,
                                                meta_data_train,
                                                group1_name,
                                                group2_name,
                                                NGENES){
    
    ## filter features
    count_matrix_train_filt <- feature_filtering(count_matrix_train)
    count_matrix_train_filt[is.na(count_matrix_train_filt)] <- 0
    
#     count_matrix_train_filt <- count_matrix_train_filt %>% edgeR::cpm() %>% t() %>% data.frame()
#     count_matrix_test_filt <- count_matrix_train_filt[rownames(count_matrix_train_filt),] %>% edgeR::cpm() %>% t() %>% data.frame()
    
#      return(list("count_matrix_train_sig" = count_matrix_train_filt, "count_matrix_test_sig" = count_matrix_train_filt))

    ## Run DESeq2 to get counts
    sig_genes <- suppressMessages(DESEQ2_feature_selection(
                            count_matrix_train_filt,
                             meta_data_train,
                             group1_name,
                             group2_name,
                            NGENES
    ))

    count_matrix_train_sig <- count_matrix_train_filt[sig_genes,] %>% edgeR::cpm() %>% t() %>% data.frame()
    count_matrix_test_sig <- count_matrix_test[sig_genes,] %>% edgeR::cpm() %>% t() %>% data.frame()    
    
    return(list("count_matrix_train_sig" = count_matrix_train_sig, "count_matrix_test_sig" = count_matrix_test_sig))
    
    
    
}

##------------------------------------
## STANDARD FEATURE FILTERING
##------------------------------------

# Feature filtering used with every feature selection algorithm

feature_filtering <- function(counts,
                              CORR_THRESH = 0.75, 
                              LOW_THRESH = 20
                             ){
    
    #####
    ## input: 
    # - counts = count matrix, normalized cpm (Columns are genes, rows are samples)
    # - CORR_THRESH = max corerlation to not remove a gene
    # - LOW_THRESH = minimum sum of all counts for that gene to be removed
    ## output:
    # - counts.cor.var.cnts = filtered count matrix
    #####
    
    #suppressMessages(library(caret))
    
    ## Convert to CPM for filtering
    counts_cpm <- counts %>% edgeR::cpm() %>% t() %>% data.frame()
    
    ## Select genes to remove based on counts
    nzv <- nearZeroVar(counts_cpm)                                          # Variance
    low.cnts <- c(1:ncol(counts_cpm))[colSums(counts_cpm) < LOW_THRESH]     # Low counts

    remove <- c(nzv,low.cnts)
    counts_cpm_sub <- counts_cpm[,-remove]
    
    ## Select genes to remove based on correlation
#     high.cor <- findCorrelation(cor(counts_cpm_sub,method = "pearson"),     # Correlation
#                                 cutoff = CORR_THRESH)   

#     counts_cpm_sub <- counts_cpm_sub[,-high.cor]
    
	thresh = dim(counts)[2]*0.75*0.5
	cut = counts %>% edgeR::cpm() %>% t() %>% colSums() %>% data.frame() %>% filter(. < thresh) %>% rownames()    
    counts_sub <- counts[colnames(counts_cpm_sub),]
    counts = counts[!(rownames(counts) %in% cut),]

    return(counts_sub)
    
}


##------------------------------------
## DESEQ2 FEATURE SELECTION
##------------------------------------

# Feature filtering used with every feature selection algorithm

DESEQ2_feature_selection <- function(counts,
                                     meta_data,
                                     group1_name,
                                     group2_name,
                                     NGENES,
                                     BASEMEAN = 0,
                                     gene_biotype_key_path = "/workdir/cfrna/references/human/hg38/gencode.biotype.name.key.tsv"
                             ){
    
    #####
    ## input: 
    # - counts = count matrix, normalized cpm (Columns are genes, rows are samples)
    # - CORR_THRESH = max corerlation to not remove a gene
    # - LOW_THRESH = minimum sum of all counts for that gene to be removed
    ## output:
    # - counts.cor.var.cnts = filtered count matrix
    #####
    
    suppressMessages(library(DESeq2))
    suppressMessages(library(data.table))
    suppressMessages(library(edgeR))
    suppressMessages(library(limma))
    ##------------------------------------
    # Contstruct DESeq Data Set
    dds <- DESeqDataSetFromMatrix(round(counts),
                                    colData = meta_data,
                                    design = ~ group + cohort +0)

    ##------------------------------------
    # Add Gene metadata
    annotation = fread(file=gene_biotype_key_path)
    annotation <- annotation[match(gsub("\\_.*","",rownames(dds)), annotation$gene_id),]
    all(rownames(dds) == annotation$ftcount_id)
    mcols(dds) <- cbind(mcols(dds), annotation)


    ##------------------------------------
    # Re-factor
    dds$group <- factor(dds$group, levels = c(group1_name, group2_name))

    ##------------------------------------
    # DAA
    dds <- DESeq(dds)
    vsd = vst(dds, blind=FALSE)
    mat = assay(vsd)
    mat = limma::removeBatchEffect(mat, vsd$cohort)
    assay(vsd) = mat
    counts_batch_corrected = assay(vsd)

    d0 = DGEList(counts_batch_corrected)
    d0 = calcNormFactors(d0)
    mm = model.matrix(~0+group+cohort, meta_data)
    y = voom(d0, mm, plot=F)
    fit = lmFit(y,mm)
    contr = makeContrasts(paste0("group",group1_name ,"- group",group2_name), levels=colnames(coef(fit)))
    tmp = contrasts.fit(fit, contr)
    tmp = eBayes(tmp)

    res = topTable(tmp, sort.by="P", n=NGENES)
    SIG_GENES = rownames(res)
    ##------------------------------------
    # Get results
    #res <- results(dds,alpha=0.01, contrast = c("group",group1_name, group2_name)) %>% data.frame() 

    #res$gene_name <- mcols(dds)$gene_name
    #res$gene_type <- mcols(dds)$gene_type

    #SIG_GENES <- res %>%
    #                filter(padj < 0.1) %>%
    #                filter(gene_type == "protein_coding") %>%
    #                filter(baseMean > BASEMEAN) %>%
    #                arrange(desc(abs(log2FoldChange))) %>%
#   #                  arrange(padj) %>%
    #                rownames() %>%
    #                head(NGENES)

    ##------------------------------------
    # Subset and normalize count matrices
    counts_sub <- counts[SIG_GENES,] 

    return(SIG_GENES)
}
