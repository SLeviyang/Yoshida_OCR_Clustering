library(dplyr)
library(plyr)
library(Biostrings)
library(SummarizedExperiment)

get_Yoshida_motifs <- function()
{
  load("../../data_input/Yoshida_TF_PWM.rda")
  return (Yoshida_motifs)
}



# Returns the number of motif occurences over a collection of sequences
motif_counts <- function(s, motifs)
{
  cc <- motifmatchr::matchMotifs(motifs, s, out="matches", p.cutoff=.01/length(motifs)) 
  m <-  SummarizedExperiment::assay(cc) %>% as.matrix
  
  return (colSums(m))
}

# Returns the number of motif occurrences over a collection of sequences after permutation
motif_counts_randomized <- function(s, motifs, N=4)
{
  sm <- sapply(strsplit(s, split=""), function(cs) cs) %>% t
  nc <- ncol(sm)
  mc <- replicate(N, apply(sm, 1, function(r) r[sample.int(nc)] %>% paste(collapse="")) %>% 
                     motif_counts(., motifs)) %>% t
           
  
  return (mc)
}

rdna <- readDNAStringSet("response.fasta") %>% as.matrix
cdna <- readDNAStringSet("control.fasta") %>% as.matrix
w <- ncol(rdna)

rseq <- apply(rdna, 1, paste, collapse="")
rseq_rand <- apply(rdna, 1, function(r) r[sample.int(w)] %>% paste(collapse=""))


tb <- rbind(data.frame(sequence=rseq,
                       cluster=1, stringsAsFactors = F),
            data.frame(sequence=rseq_rand,
                       cluster=2, stringsAsFactors = F))
#' 
#' 
#' @param sequences a data.frame with fields sequence (characters) and
#' cluster (int)
#' @param motifs a TFBSTools PWMatrixList 
motifCluster <- function(sequences, clusters, motifs)
{
  return (list(sequences=sequences,
               clusters=clusters,
               motifs=motifs))
}

motifcounts.motifCluster <- function(mc, response_cluster)
{
  rseq <- mc$sequences[mc$clusters == response_cluster]
  cseq <- mc$sequences[mc$clusters != response_cluster]
  
  rmatch <- motifmatchr::matchMotifs(mc$motifs, rseq, out="matches", p.cutoff=1E-6)
  cmatch <- motifmatchr::matchMotifs(mc$motifs, cseq, out="matches", p.cutoff=1E-6)
  
  rm <- SummarizedExperiment::assay(rmatch) %>% as.matrix 
  cm <- SummarizedExperiment::assay(cmatch) %>% as.matrix 
  
  Nr <- nrow(rm); success_r <- colSums(rm)
  Nc <- nrow(cm); success_c <- colSums(cm)
 
  
  browser()
}