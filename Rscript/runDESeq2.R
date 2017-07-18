#!/usr/bin/env Rscript
# -- Set up the command line arguments
###########
suppressPackageStartupMessages(library("argparser"))

p <- arg_parser(description = "Submodule 3c - Run differential expression.")
p <- add_argument(p, arg = "--directory", default="NULL",
                  help="the directory where the transcript.counts, lib_data and output files are stored")
p <- add_argument(p, arg = "--counts", default="processed.counts",
                  help="the name of the transcript.counts file")
p <- add_argument(p, arg = "--lib_info", default="lib_info.csv",
                  help="the library information dataset.")
p <- add_argument(p, arg = "--baseFormula", default="y ~ 1",
                  help="The full formula with all factors to test")
p <- add_argument(p, arg = "--testName", default="test",
                  help="The name to append to columns of stats")
p <- add_argument(p, arg = "--redFormula",  default="y ~ 0",
                  help="the formula that lacks the term to test")
p <- add_argument(p, arg = "--subsetColumn",  default="NULL",
                  help="the name of the column that can be subset. If specified, so must --factor2keep")
p <- add_argument(p, arg = "--factors2keep",  default="NULL",
                  help="a vector of factor levels that will be retained. If specified, so must --subsetColumn")
p <- add_argument(p, arg = "--outputfilename",  default="de.stats",
                  help="name of stats file name")

args <- parse_args(p)

if(args$directory == "NULL") args$directory<-getwd()
if(args$factors2keep != "NULL") args$factors2keep<-strsplit(args$factors2keep, ",")[[1]]
if(any(grepl(",",c(args$baseFormula,args$redFormula, args$testName)))){
  args$baseFormula<-strsplit(args$baseFormula, ",")[[1]]
  args$redFormula<-strsplit(args$redFormula, ",")[[1]]
  args$testName<-strsplit(args$testName, ",")[[1]]
  if(length(unique(length(args$baseFormula),
                   length(args$redFormula),
                   length(args$testName)))!=1)
    stop("if multiple formulae are supplied ...
         --testName --redFormula --baseFormula must all be the same length\n")
}

###########

# -- Set up the Rscript environment
###########
suppressPackageStartupMessages(library("deTools"))
suppressPackageStartupMessages(library("DESeq2"))
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

# -- Run differential expression analysis
###########
cat("analysis includes", nrow(counts), "genes and",ncol(counts),"libraries...\n")
stats<-pipeDESeq2(counts=counts, info=info,
                  formula = args$baseFormula,
                  reduced = args$redFormula,
                  testNames = args$testName)

# -- Output results
write.table(stats$LRT_results,
            file = paste0(args$directory,"/",args$outputfilename), sep = "\t", row.names = F)
###########
