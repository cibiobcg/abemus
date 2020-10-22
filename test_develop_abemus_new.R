library( devtools )
devtools::install_github("cibiobcg/abemus", build_vignettes = T,ref = 'develop')

library(data.table)
library(parallel)
library(Matrix)

outdir <-  "/BCGLAB/ncasiraghi/abemus_test_develop_branch"
sample.info.file <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/sample_info.txt"
targetbed <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/TST170_DNA_target.nochr.merged.nochromY.bed"
pacbamfolder_bychrom <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/bychrom/"

sif <- import_sif(sample.info.file)

bed <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)

getInfo <- function(chrom,sif,pacbamfolder_bychrom,select){
  germlineset <- get_germlineset(sifgerm = sif$df_ctrl,
                                 pacbamfolder_bychrom = pacbamfolder_bychrom,
                                 chrom = chrom)
  df <- fread(germlineset[1],sep = '\t',stringsAsFactors = FALSE,skip = 1,select = select,verbose = FALSE,data.table=FALSE, nThread = 4)
  return(as.vector(unlist(df)))
}

pileup2mat <- function(chrom,sif,pacbamfolder_bychrom,select,sparse=FALSE){
  germlineset <- get_germlineset(sifgerm = sif$df_ctrl,
                                 pacbamfolder_bychrom = pacbamfolder_bychrom,
                                 chrom = chrom)
  df <- lapply(germlineset,fread,sep = '\t',stringsAsFactors = FALSE,skip = 1,select = select,verbose = FALSE,data.table=FALSE, nThread = 4)
  if(sparse){
    return(Matrix(as.matrix(do.call(cbind,df)), sparse = TRUE))
  } else {
    return(as.matrix(do.call(cbind,df)))
  }
}

ref <- lapply(bed$chromosomes,getInfo,sif,pacbamfolder_bychrom,select=3)
pos <- lapply(bed$chromosomes,getInfo,sif,pacbamfolder_bychrom,select=2)

mat_vaf <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=8,sparse=TRUE)
mat_cov <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=9)

dict <- 4:7
names(dict) <- c('A','C','G','T')
mat <- list()
for(base in names(dict)){
  mat[[base]] <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=as.integer(dict[base]),sparse=TRUE)
}

# compute pbem

# select correct rows!

num <- rowSums((mat[['A']][[6]] + mat[['C']][[6]] + mat[['G']][[6]]))
den <- rowSums(mat_cov[[6]])

pbem <- num/den

hist(pbem,100)

# compute af threshold








outdir.bperr.name = "BaseErrorModel"
coverage_binning = 50
af_max_to_compute_thresholds = 0.2
coverage_min_to_compute_thresholds = 10
af_max_to_compute_pbem = 0.2
coverage_min_to_compute_pbem = 10
n_pos_af_th = 0.2
mc.cores = 1
step = 5000



sif <- read.delim(file = sample.info.file,header = F)

sif$V5 <- file.path('/BCGLAB/ncasiraghi/abemus_test_develop_branch/bychrom/')

bed <- read.delim(file = targetbed,header = F)

head(bed)

bed <- bed[which(bed$V1 != 'Y'),]
write.table(bed,file = "/BCGLAB/ncasiraghi/abemus_test_develop_branch/TST170_DNA_target.nochr.merged.nochromY.bed",quote = F,sep = '\t',col.names = F,row.names = F)
