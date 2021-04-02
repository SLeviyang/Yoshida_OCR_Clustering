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

compute_chromVar_deviations <- function(args)
{
  count_file <- args[1]
  motif_file <- args[2]
  bed_file <- args[3]
  out_file <- args[4]
  
  counts <- data.table::fread(count_file, sep=",", header=T) %>% as.matrix()
  m <- data.table::fread(motif_file, sep=",", header=T) %>% as.matrix()
  b <- data.table::fread(bed_file, sep="\t", header=F) %>%
       setNames(c("chr", "chrStart", "chrEnd"))
  
  if (nrow(b) != nrow(m))
    stop("m and b should have the same dimension!")
  
  ir <- IRanges(start=b$chrStart, end=b$chrEnd)
  gr <- GRanges(seqnames=b$chr, ranges = ir, strand = "+")
  
  counts <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                 rowRanges = gr)
  counts_GC <- addGCBias(counts, genome = BSgenome.Mmusculus.UCSC.mm10)
  dev <- computeDeviations(object = counts_GC, annotations = m)
  z <- assays(dev)$deviations
  colnames(z) <- colnames(counts)
  d <- data.frame(motif=colnames(m), stringsAsFactors = F) %>%
       cbind(z)
  
  write.csv(d, out_file, row.names=F)
}

compute_chromVar_deviations(args)



