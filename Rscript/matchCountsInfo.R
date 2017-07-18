#!/usr/bin/env Rscript
# -- Set up the command line arguments
###########
suppressPackageStartupMessages(library("argparser"))

p <- arg_parser(description = "Submodule 3b - match counts data with library information data. Takes any counts-formatted dataset and matches the columns names with the library info data.")
p <- add_argument(p, arg = "--directory", default="NULL",
                  help="the directory where the counts, lib_data.txt and output files are stored")
p <- add_argument(p, arg = "--lib_data", default = "lib_info.csv",
                  help="the name of the lib_data.txt file")
p <- add_argument(p, arg = "--counts", default = "processed.counts",
                  help="the name of the .counts file")
p <- add_argument(p, arg = "--idcolumn", default = "library",
                  help="the name of the .counts file")

args <- parse_args(p)

if(args$directory == "NULL") args$directory<-paste0(getwd(),"/")
idcol<-args$idcolumn
###########

# -- Read in the counts data
###########
counts<-read.delim(paste0(args$directory, args$counts), header=T)
info<-read.csv(paste0(args$directory, "/",args$lib_data), header=T, stringsAsFactors = F)

# -- Find the matching columns in the library data with
if(idcol %in% colnames(info) & sum(info[,idcol] %in% colnames(counts))>1){
  colnames(info)[which(colnames(info)=="idcol")]<-"library"
  ids<-intersect(colnames(counts),info$library)
}else{
  noverl<-sapply(info, function(x) sum(x %in% colnames(counts)))
  colnames(info)[which.max(noverl)]<-"library"
  ids<-intersect(colnames(counts),info$library)
}
rownames(info)<-info$library
count.out<-counts[,ids]
info<-info[ids,]
write.csv(info, file = paste0(args$directory, "/lib_info.csv"), row.names = F)
write.table(count.out, file = paste0(args$directory, "/processed.counts"), sep = "\t")
###########
