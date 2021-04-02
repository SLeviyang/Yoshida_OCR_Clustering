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

INPUT_DATA_DIR = "../../data_input/"
DATA_DIR <- "../../data_output/"

get_Yoshida_motifs <- function()
{
  data("mouse_pwms_v2")
  motifs <- mouse_pwms_v2
  names(motifs) <- sapply(motifs, name) %>% toupper
  
  Yoshida_TF_file <- paste(INPUT_DATA_DIR,
                           "Yoshida_Table_S5.csv",   
                           sep="")
  TF <- read.csv(Yoshida_TF_file, header=T, stringsAsFactors = F) %>% dplyr::pull(TF) %>% toupper
  motifs <- motifs[is.element(names(motifs), TF)]
  
  return (motifs)
}

load_workflow <- function(binary_matrix_source,
                                        edge_FDR,
                                        idr_FDR = NULL,
                                        count_cutoff = NULL)
{
  if (binary_matrix_source == "count")
    tag <- "Yoshida_count_matrix_cutoff_" + str(count_cutoff) 
  else {
      tag <-  paste("Yoshida_idr_matrix_idr_", idr_FDR, sep="")
      tag = paste(tag, "_edge_FDR_", edge_FDR, sep="")
  }
      
  print(c("accessing workflow", tag))
      
  output_dir = paste(DATA_DIR, tag, "/", sep="")
  # DATA_DIR defined in configuration.R
  cluster_file <- paste(output_dir, "clusters.csv", sep="")
  sequence_file <- paste(output_dir, "sequences.csv", sep="")
  matrix_file <- paste(output_dir, "matrix.csv", sep="")
  genomics_file <- paste(output_dir, "genomics.csv", sep="")
  bed_file <- paste(output_dir, "loci.bed", sep="")
  
  pc <- list(clusters=read.csv(cluster_file, header=T, stringsAsFactors = F),
             sequences=read.csv(sequence_file, header=F, stringsAsFactors = F) %>%
                        dplyr::pull("V1"),
             m=read.csv(matrix_file, header=T, stringsAsFactors = F),
             g=read.csv(genomics_file, header=T, stringsAsFactors = F),
             b=read.table(bed_file, sep="\t", header=F) %>%
               setNames(c("chr", "chrStart", "chrEnd")),
             output_dir=output_dir)
  
  return (pc)
}

create_motif_matrix <- function(w)
{
  motifs <- get_Yoshida_motifs()
  outfile <- paste(w$output_dir, "motif.csv", sep="")
 
  m <- w$m %>% as.matrix
  seqs <- DNAStringSet(w$sequences)
  motif_ix <- matchMotifs(motifs, seqs, bg="subject", out="matches",
                          p.cutoff=5E-6)
  
  motif_matches <- assays(motif_ix)$motifMatches %>% as.matrix
  data.table::fwrite(motif_matches, outfile, sep=",", row.names=F)
  
  return (NULL)
}

# create matrix that maps peaks to peaks with similar GC and counts/accessibility
create_background_matrix <- function(w)
{
  m <- w$m %>% as.matrix
  nonzero_rows <- which(rowSums(m) > 0)
  
  m <- m[nonzero_rows,]
  b <- w$b[nonzero_rows,]
  ir <- IRanges(start=b$chrStart, end=b$chrEnd)
  gr <- GRanges(seqnames=b$chr, ranges = ir, strand = "+")
  
  counts <- SummarizedExperiment(assays=SimpleList(counts=m),
                                 colData=data.frame(cell_type=colnames(m),
                                                    stringsAsFactors = F),
                                 rowRanges = gr)
  counts_GC <- addGCBias(counts, 
                         genome = BSgenome.Mmusculus.UCSC.mm10)
  bg <- getBackgroundPeaks(object = counts_GC) 
  colnames(bg) <- paste("map", 1:ncol(bg), sep="")
  d <- data.frame(row=nonzero_rows) %>% cbind(bg)
  
  outfile <- paste(w$output_dir, "peak_chromVar_maps.csv", sep="")
  data.table::fwrite(d, outfile, sep=",", row.names=F)
  
  return (NULL)
}



