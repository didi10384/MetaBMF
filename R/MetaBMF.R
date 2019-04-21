###Load Dependent Packages
library('seqinr')
library('getopt')
library('MASS')

###Get Options
myopt <- matrix(c('help','h',0,"logical",
                 'work_dir','w',1,"character",
                 'metamat_path','m',1,"character",
                 'num_cup_cores','n',2,"integer",
                 'bin_step','e',1,"integer",
                 'bin_min','l',1,"integer",
                 'bin_max','u',1,"integer",
                 'ctg_len_trim','r',2,'integer'),byrow = TRUE,ncol = 4)

opt <- getopt(myopt)

if( !is.null(opt$help)) {
  cat(getopt(myopt, usage=TRUE));
  q(status=1);
}

if ( is.null(opt$num_cup_cores ) ) { opt$num_cup_cores = 1 }
#if ( is.null(opt$bin_step ) ) { opt$bin_step = 10 }
#if ( is.null(opt$bin_min ) ) { opt$bin_min = 10 }
#if ( is.null(opt$bin_max ) ) { opt$bin_max = 1000 }
if ( is.null(opt$ctg_len_trim ) ) { opt$ctg_len_trim = 1000 }
#if ( is.null(opt$plot_CK ) ) { opt$plot_bic = "F" }


work_dir <- opt$work_dir
metamat_path<- opt$metamat_path

bin_step <- opt$bin_step
bin_min <- opt$bin_min
bin_max <- opt$bin_max

num_cpu <- opt$num_cup_cores
ctg_len_trim <- opt$ctg_len_trim

out_dir <- paste(work_dir,"/output",sep="")

source(paste(metamat_path,"/R/seeds.R",sep=""))
source(paste(metamat_path,"/R/silhouette.R",sep=""))
source(paste(metamat_path,"/R/cluster.R",sep=""))
source(paste(metamat_path,"/R/abundance.R",sep=""))

###Load the Input Dataset
rcmm_file <- paste(work_dir,"/output/count-map.tsv",sep="")
ctg_file <- paste(work_dir,"/contigs/Contigs.fasta",sep="")
dmat <- as.matrix(read.table(rcmm_file,header=T,check.names = FALSE))
fadat <- read.fasta(file = ctg_file)

ctg_name <- unlist(lapply(getAnnot(fadat), function(x){strsplit(x," ")[[1]][1]}))

lvec <- as.numeric(unlist(lapply(fadat,length)))

nvec <- numeric(ncol(dmat))

reads_sum <- read.table(paste(work_dir,"/reads_info.txt",sep=""), sep=",", header=T)
nvec <- reads_sum[,3]
sample_name <- reads_sum[,1]

if(length(which(nvec>0.1))!=length(sample_name)){
  cat("Error: sample_name are not matched.");
  cat(sample_name)
  q(status=1);
}

###Preprocessing the data
dmat_rsum <- apply(dmat,1,sum)
ind_in <- which(lvec>=ctg_len_trim&dmat_rsum>=200)
dmat_in <- dmat[ind_in,]

###Select Seed Contigs
Seeds <- seeds(dmat_in,bin_max)

###Select Num of Species using Silhouette Statistic
pmat <- Seeds$pmat
dmat_n <- Seeds$dmat_n

ks.i <- NULL
for(i in seq(bin_min,bin_max,bin_step))
{
  ks.i <- rbind(ks.i,c(i,silhouette(dmat_n,pmat,i)))
}

K <- ks.i[which.max(ks.i[,2]),1]
cat("The Optimal Number of Species is",K,"\n")

###Output the Cluster Label
segs <- cluster(dmat_n,pmat,K)$clus

write.table(cbind(ctg_name[ind_in], segs), file=paste(work_dir,'/output/segs.txt',sep=""), col.names=F, row.names=F)
cat("Succesffully Output the Binning Results.\n")

###Output the Relative Abundance Matrix
pmat1 <- Seeds$pmat1[1:K,]
sample_name = colnames(dmat)
xmat <- t(abundance(pmat1,dmat_in,segs,lvec[ind_in],nvec))
colnames(xmat) <- sample_name

write.table(xmat, file=paste(work_dir,'/output/relative_abundance.txt',sep=""), col.names=T, row.names=F)
cat("Succesffully Output the Relative Abundance Matrix.\n")
