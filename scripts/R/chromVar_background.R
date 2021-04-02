library(dplyr)
library(plyr)
library(Biostrings)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVARmotifs)
library(TFBSTools)

args = commandArgs(trailingOnly=TRUE)
print(args)

create_background_matrix <- function(args)
{
  X_file <- args[1]
  bed_file <- args[2]
  out_file <- args[3]
  
  m <- data.table::fread(X_file, sep=",", header=F) %>% as.matrix()
  b <- data.table::fread(bed_file, sep="\t", header=F) %>%
       setNames(c("chr", "chrStart", "chrEnd"))
  
  if (nrow(b) != nrow(m))
    stop("m and b should have the same dimension!")
  
  ir <- IRanges(start=b$chrStart, end=b$chrEnd)
  gr <- GRanges(seqnames=b$chr, ranges = ir, strand = "+")
  
  counts <- SummarizedExperiment(assays=SimpleList(counts=m),
                                 rowRanges = gr)
  counts_GC <- addGCBias(counts, genome = BSgenome.Mmusculus.UCSC.mm10)
  bg <- getBackgroundPeaks(object = counts_GC) 
  
  data.table::fwrite(bg, out_file, sep=",", row.names=F, col.names = F)
}

create_background_matrix(args)
