#!/usr/bin/env Rscript
# -- Set up the command line arguments
###########
suppressPackageStartupMessages(library("argparser"))

p <- arg_parser(description = "Submodule 3d - Calculate voom-normalized counts.")
p <- add_argument(p, arg = "--directory", default="NULL",
                  help="the directory where the transcript.counts, lib_data and output files are stored")
p <- add_argument(p, arg = "--counts", default="processed.counts",
                  help="the name of the transcript.counts file")
p <- add_argument(p, arg = "--lib_info", default="lib_info.csv",
                  help="the library information dataset.")
p <- add_argument(p, arg = "--subsetColumn",  default="NULL",
                  help="the name of the column that can be subset. If specified, so must --factor2keep")
p <- add_argument(p, arg = "--factors2keep",  default="NULL",
                  help="a vector of factor levels that will be retained. If specified, so must --subsetColumn")
p <- add_argument(p, arg = "--outputfilename",  default="voom.counts",
                  help="name of voom counts file name")

args <- parse_args(p)

if(args$directory == "NULL") args$directory<-getwd()
if(args$factors2keep != "NULL") args$factors2keep<-strsplit(args$factors2keep, ",")[[1]]

###########

# -- Set up the Rscript environment
###########
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
# -- Read in the counts data
###########
counts<-read.delim(paste0(args$directory, "/",args$counts), header=T)
info<-read.csv(paste0(args$directory, "/",args$lib_info), header=T, stringsAsFactors = F)
if(!identical(info$library, colnames(counts)))
  stop("library and counts data do not match - run matchCountsInfo.R first\n")
###########

# -- Cull to subset (if necessary)
###########
if(all(args$subsetColumn != "NULL") & all(args$factors2keep != "NULL")){
  if(!args$subsetColumn %in% colnames(info))
    stop("subsetColumn is not a column in lib_info")
  if(all(args$factors2keep == "NULL") | all(args$subsetColumn == "NULL"))
    stop("if subsetColumn is provided, so must factors2keep")
  wh<-which(info[,args$subsetColumn] %in% args$factors2keep)
  counts<-counts[,wh]
  info<-info[wh,]
}

# -- Calc voom
###########
v<-voomWithQualityWeights(counts = DGEList(counts = counts))

# -- Output results
write.table(v$E,
            file = paste0(args$directory,"/",args$outputfilename),
            sep = "\t", row.names = F)
###########
