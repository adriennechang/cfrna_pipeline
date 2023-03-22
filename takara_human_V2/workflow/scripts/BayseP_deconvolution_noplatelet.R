suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(BayesPrism))

###----------------------------------------------
# Inputs
args<-commandArgs(TRUE)

cfrna_counts <- args[1]
ref_counts <- args[2]
ref_meta <- args[3]
biotype_key <- args[4]
output_file <- args[5]
output_file_pc <- args[6]
output_file_pc_sig <- args[7]
cores <- as.numeric(args[8])
sample_id <- args[9]

###----------------------------------------------
# Read in data
cat("-----> READING IN DATA \n")

rm_platelet = c('ENSG00000162511.8','ENSG00000117400.18','ENSG00000162366.8','ENSG00000163220.11','ENSG00000242252.2',
                'ENSG00000066294.15','ENSG00000158769.18','ENSG00000143226.15','ENSG00000198734.12','ENSG00000174175.17','ENSG00000073756.13','ENSG00000150681.10','ENSG00000117335.20','ENSG00000116977.19','ENSG00000205639.12','ENSG00000228474.6','ENSG00000115956.10','ENSG00000115008.6','ENSG00000121966.7','ENSG00000115935.18','ENSG00000168497.5','ENSG00000144476.6','ENSG00000183813.7','ENSG00000183625.16','ENSG00000163823.4','ENSG00000233276.8','ENSG00000163932.15','ENSG00000169704.5','ENSG00000169313.10','ENSG00000178732.5','ENSG00000174125.8','ENSG00000174130.12','ENSG00000124406.16','ENSG00000163737.4','ENSG00000163736.4','ENSG00000138722.10','ENSG00000155016.18','ENSG00000071205.12','ENSG00000129116.19','ENSG00000153071.15','ENSG00000164171.11','ENSG00000181104.7','ENSG00000145824.13','ENSG00000113140.11','ENSG00000124491.16','ENSG00000010704.19','ENSG00000206503.13','ENSG00000204420.10','ENSG00000204310.13','ENSG00000112062.11','ENSG00000161911.11','ENSG00000112195.9','ENSG00000171611.10','ENSG00000112715.25','ENSG00000156535.15','ENSG00000198478.8','ENSG00000135604.10','ENSG00000135218.19','ENSG00000127920.6','ENSG00000106366.9','ENSG00000005249.13','ENSG00000091137.14','ENSG00000104267.10','ENSG00000154188.10','ENSG00000137033.12','ENSG00000136869.16','ENSG00000095303.17','ENSG00000266412.6','ENSG00000122862.5','ENSG00000107438.9','ENSG00000142082.15','ENSG00000166333.14','ENSG00000149781.13','ENSG00000068831.19','ENSG00000137507.11','ENSG00000154146.13','ENSG00000149564.12','ENSG00000110799.14','ENSG00000010278.15','ENSG00000111321.11','ENSG00000059804.16','ENSG00000110848.8','ENSG00000165682.14','ENSG00000111348.9','ENSG00000110934.13','ENSG00000127314.18','ENSG00000110876.10','ENSG00000111275.13','ENSG00000073060.17','ENSG00000132965.10','ENSG00000125257.16','ENSG00000057593.14','ENSG00000197930.13','ENSG00000137801.11','ENSG00000137845.15','ENSG00000255346.10','ENSG00000167972.14','ENSG00000197471.12','ENSG00000167207.14','ENSG00000074370.18','ENSG00000185245.8','ENSG00000108839.12','ENSG00000271503.6','ENSG00000005961.19','ENSG00000259207.9','ENSG00000261371.6','ENSG00000150637.9','ENSG00000125735.11','ENSG00000129355.7','ENSG00000123146.20','ENSG00000088826.18','ENSG00000101335.10','ENSG00000064601.21','ENSG00000101017.14','ENSG00000124126.14','ENSG00000101162.4','ENSG00000101194.18','ENSG00000125534.10','ENSG00000156265.16','ENSG00000160255.18','ENSG00000203618.6','ENSG00000100351.17','ENSG00000205755.12','ENSG00000212747.5','ENSG00000022267.19','ENSG00000102245.8','ENSG00000205755.12_PAR_Y')

bk.dat = read.delim(cfrna_counts,row.names=c(1),skip=1)
bk.dat <- bk.dat[,7,drop= FALSE]
colnames(bk.dat) <- sample_id

bk.dat <- bk.dat %>% t() %>% data.frame()


sc.dat <- as.matrix(fread(ref_counts),rownames=1)

sc.meta <- fread(ref_meta)
cell.type.labels = sc.meta$cell_type_collapsed

###----------------------------------------------
# Clean up reference matrix

gene_group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY")

cat("-----> CLEANING UP REFERENCE \n")
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                    species="hs", 
                                    gene.group=gene_group ,
                                    exp.cells=5)
     

###----------------------------------------------
# ALL GENES
###----------------------------------------------

###----------------------------------------------
# Create prism                         
cat("-----> CREATE PRISM \n")                 
myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)
          
###----------------------------------------------
# Run deconvolution                          
bp.res <- run.prism(prism = myPrism, n.cores=cores)
theta <- get.fraction (bp=bp.res,
            which.theta="final",
            state.or.type="type")

###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")     
output_theta = data.frame(theta)
output_theta$sample_id <- rownames(output_theta)
output_theta = output_theta %>% select(sample_id, everything())
write.table(output_theta,output_file,sep="\t",quote=FALSE,row.names=FALSE)



###----------------------------------------------
# PROTEIN CODING ONLY 
###----------------------------------------------

###----------------------------------------------
# Extract only protein coding genes
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")

###----------------------------------------------
# Create prism                         
cat("-----> CREATE PRISM \n")                 
myPrism.pc <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)
          
###----------------------------------------------
# Run deconvolution                          
bp.pc.res <- run.prism(prism = myPrism.pc, n.cores=cores)
theta.pc <- get.fraction (bp=bp.pc.res,
            which.theta="final",
            state.or.type="type")

###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")    
output_theta.pc = data.frame(theta.pc)
output_theta.pc$sample_id <- rownames(output_theta.pc)  
output_theta.pc = output_theta.pc %>% select(sample_id, everything())              
write.table(output_theta.pc,output_file_pc,sep="\t",quote=FALSE,row.names=FALSE)
    




###----------------------------------------------
# MARKER GENES & PROTEIN CODING ONLY 
###----------------------------------------------

###----------------------------------------------
# Perform differential expression
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat.filtered.pc[,colSums(sc.dat.filtered.pc>0)>3],# filter genes to reduce memory use
                                          cell.type.labels=cell.type.labels,
                                          cell.state.labels=cell.type.labels,
                                          psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                          cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                          n.cores=1 #number of threads
)

###----------------------------------------------
# Extract protein coding genes
sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                                  stat=diff.exp.stat,
                                                  pval.max=0.01,
                                                  lfc.min=0.1)

###----------------------------------------------
# Create prism
myPrism.pc.sig <- new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1
)

###----------------------------------------------
# Run deconvolution    
bp.pc.sig.res <- run.prism(prism = myPrism.pc.sig, n.cores=50)

theta.pc.sig <- get.fraction (bp=bp.pc.sig.res,
            which.theta="final",
            state.or.type="type")

###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")    
output_theta.pc.sig = data.frame(theta.pc.sig)
output_theta.pc.sig$sample_id <- rownames(output_theta.pc.sig)  
output_theta.pc.sig = output_theta.pc.sig %>% select(sample_id, everything())              
write.table(output_theta.pc.sig,output_file_pc_sig,sep="\t",quote=FALSE,row.names=FALSE)