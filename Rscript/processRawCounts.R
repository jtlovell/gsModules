#!/usr/bin/env Rscript
# -- Set up the command line arguments
###########
suppressPackageStartupMessages(library("argparser"))

p <- arg_parser(description = "Submodule 3b - process raw RNA-seq counts data. Takes a raw transcript abundance (tab delimited) matrix where rows are genes and columns are libraries. Returns a matrix where libraries that are outliers (too large, too small) and genes that are not highly expressed are dropped.")
p <- add_argument(p, arg = "--directory", default="NULL",
                  help="the directory where the transcript.counts and output files are stored")
p <- add_argument(p, arg = "--counts", default="transcript.counts",
                  help="the name of the transcript.counts file")
p <- add_argument(p, arg = "--libsizeoutlier", default=0.01,
                    help="the significance level at which library size outliers are dropped")
p <- add_argument(p, arg = "--threshold", default=5,
                    help="The minimum mean counts / gene to retain (converted to cpm)")
p <- add_argument(p, arg = "--output",  default="processed.counts",
                    help="the name of the processed counts file")

args <- parse_args(p)

if(args$directory == "NULL") args$directory<-getwd()
###########

# -- Set up the Rscript environment
###########
suppressPackageStartupMessages(library("outliers"))
suppressPackageStartupMessages(library("edgeR"))

outl.grubbs<-function(x, plotit=T, thresh = 0.01){
  sig<-T
  while(sig){
    ot<-outliers::grubbs.test(x)
    otl.num<-as.numeric(strsplit(ot[[2]], " ")[[1]][3])
    index<-which(x == otl.num & !is.na(x))[1]
    if(ot$p.value<=thresh){
      x[index]<-NA
    }else{
      sig<-F
    }
  }
  return(x)
}
###########

# -- Read in the counts data
###########
cat("reading in counts dataset with ...")
counts<-read.delim(paste0(args$directory,"/",args$counts), header=T)
cat(nrow(counts),"genes and", ncol(counts), "libraries\n")
counts<-data.matrix(counts)

###########

# -- Drop very small or very large libraries
###########
ol<-outl.grubbs(colSums(counts), thresh = args$l)
todrop<-which(is.na(ol))
cat("dropping libs:", paste(colnames(counts)[todrop], collapse = ", "),"\n")
counts<-counts[,-todrop]
###########

# -- Determine counts / million
###########
countsPM<-cpm(counts)
geneCountsPM<-rowMeans(countsPM)
geneCounts<-rowMeans(counts)
###########

# -- Determine threshold based on cpm - raw count relationship
###########
mod<-lm(geneCountsPM ~ geneCounts)
slope<-summary(mod)$coefficients[2,"Estimate"]
cpmThresh <- slope*args$t
###########

# -- Drop genes with mean cpm < threshold and write files
###########
counts<-counts[geneCountsPM>=cpmThresh,]
cat("retaining", nrow(counts),"genes and",ncol(counts),"libraries\n")
write.table(counts, file = paste0(args$directory, "processed.counts"), sep = "\t")
###########
