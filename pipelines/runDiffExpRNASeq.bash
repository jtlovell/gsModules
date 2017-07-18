# sm 3a get fastq processing
# sm 3b make counts
# sm 3c prepare counts for differential expression
## retain libraries that are not significant (P <= 0.001) outliers in terms of read counts
## retain genes that have ~5 raw counts
## if the counts file is called transcript.counts, does not need to be specified
./processRawCounts.R \
  --directory /Users/John/Desktop/gsModules/data/rnaseq_counts \
  --libsizeoutlier 0.001 \
  --threshold 5

# sm 3c match counts and info data
./matchCountsInfo.R --directory /Users/John/Desktop/gsModules/data/rnaseq_counts

# sm 3d run differential expression
## for a matched counts and info dataset, parse so that only factors desired are retained
## run differential expression analysis using likelihood ratio tests between of hierarchical models
./runDESeq2.R \
  --directory  /Users/John/Desktop/gsModules/data/rnaseq_counts \
  --baseFormula ~id+Treatment,~id+Treatment,~id+Treatment+id*Treatment \
  --redFormula ~id,~Treatment,~id+Treatment \
  --testName trt,id,trt.id \
  --subsetColumn id \
  --factors2keep HAL2,FIL2

  # sm 3e calculate voom-normalized data
  ## same arguments as above, but without formula (since normalization occurs agnostic of experimental design)
./calcVoom.R \
  --directory  /Users/John/Desktop/gsModules/data/rnaseq_counts \
  --subsetColumn id \
  --factors2keep HAL2,FIL2

