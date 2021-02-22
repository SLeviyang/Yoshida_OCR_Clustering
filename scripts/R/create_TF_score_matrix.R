library(dplyr)
library(plyr)
library(Biostrings)
library(JASPAR2018)

source("utilities_HOMER.R")

INPUT_DATA_DIR = "../../data_input/"
DATA_DIR <- "../../data_output/"
HOMER_PATH <- "/Users/sivanleviyang/Dropbox/OSX_software/homer/bin"
TF_MOTIF_FILE <- paste(DATA_DIR, "TF_motifs.txt", sep="")

create_Yoshida_motif_file <- function()
{
  full_JASPAR_meme_file <- paste(INPUT_DATA_DIR,
                                 "JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt",
                                 sep="")
  homer_file <- paste("temp", round(1E10*runif(1)), ".txt", sep="")
  meme2homer.BM(full_JASPAR_meme_file, homer_file, min_nsites=0, set_score=-1E10)
  jaspar_motifs <- read_homer_motifs.BM(homer_file)
  full_motif_names <- get_motif_names.BM(jaspar_motifs) %>% toupper
  motif_names <- sapply(strsplit(full_motif_names, split="_"), "[", 2)
  names(jaspar_motifs) <- motif_names
  
  Yoshida_TF_file <- paste(INPUT_DATA_DIR,
                           "Yoshida_Table_S5.csv",   
                           sep="")
  TF <- read.csv(Yoshida_TF_file, header=T, stringsAsFactors = F) %>%
        dplyr::filter(Significant) %>%
        dplyr::pull(TF) %>%
        toupper %>% sort
  TF <- c(TF, "STAT1::STAT2")
  
  joint_motifs <- jaspar_motifs[is.element(names(jaspar_motifs), TF)]
  write_homer_motifs.BM(joint_motifs, TF_MOTIF_FILE)
                       
  
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
  if (!file.exists(TF_MOTIF_FILE))
    create_Yoshida_motif_file()
  
  pc <- list(clusters=read.csv(cluster_file, header=T, stringsAsFactors = F),
             sequences=read.csv(sequence_file, header=F, stringsAsFactors = F) %>%
                        dplyr::pull("V1"),
             motif_file=TF_MOTIF_FILE,
             output_dir=output_dir)
  
  return (pc)
}

create_all_cluster_TF_score_matrices <- function(edge_FDR,
                                                idr_FDR,
                                                max_cluster_num=22,
                                                ncore=4,
                                                create_null=F)
{
  w <- load_workflow("idr", edge_FDR, idr_FDR=idr_FDR)
  outdir <- paste(w$output_dir, "TF_scores/", sep="")
  if (!dir.exists(outdir))
    dir.create(outdir)

  for (cluster_num in 0:max_cluster_num) {
    if (create_null)
      outfile <- paste(outdir, "TF_null_score_matrix_cluster_", cluster_num, ".csv", sep="")
    else 
      outfile <- paste(outdir, "TF_score_matrix_cluster_", cluster_num, ".csv", sep="")
    
    if (file.exists(outfile)) 
      cat(outfile, "already exists, skippping...\n")
    else {
      cat("constructing score matrix for cluster", cluster_num, "\n")
      create_cluster_TF_score_matrix(w, cluster_num, outfile, ncore=4,
                                     create_null=create_null)
    }
  }
}

# Create a score matrix over al sequences in the cluster (rows) and
# all motifs (columns).  If create_null, then sequences are shuffled to
# create a null distribution of scores.  This is a hack because the null
# distribution can be computed approximately by binning possible scores
# and using dynamics programming
create_cluster_TF_score_matrix <- function(w, cluster_num, 
                                           outfile,
                                           ncore=4,
                                           create_null)
{
  rows <- dplyr::filter(w$clusters, cluster==cluster_num) %>%
          dplyr::pull(row)
  # adjust for python count
  rows <- rows + 1
  
  seqs <- w$sequences[rows]
  if (create_null) 
    seqs <- sapply(strsplit(seqs, split=""), function(ss) 
      paste(ss[sample.int(length(ss))], collapse=""))
 
  
  sf <- paste("temp_seq", round(1E10*runif(1)), ".csv", sep="")
  con <- file(sf, open="w")
  writeLines(seqs, con)
  close(con)
  
  create_motif_score_matrix(sequence_file=sf,
                            motif_file=w$motif_file,
                            outfile=outfile,
                            ncore=ncore,
                            HOMER_PATH=HOMER_PATH)
  file.remove(sf)
}


###################################################
# work functions
# load_motifs <- function(meme_file)
# {
#   homer_file <- paste("temp", round(1E10*runif(1)), ".txt", sep="")
#   meme2homer.BM(meme_file, homer_file, min_nsites=0, set_score=-1E10)
#   
#   motifs <- read_homer_motifs.BM(homer_file)
#   file.remove(homer_file)
#   
#   full_motif_names <- get_motif_names.BM(motifs) %>% toupper
#   motif_names <- sapply(strsplit(full_motif_names, split="_"), "[", 2)
#   
#   if (length(motif_names) != length(unique(motif_names)))
#     stop("motif names are not unique!")
#   names(motifs) <- motif_names
#   
#   return (motifs)
# }


#' @param sequence_file a file containing sequences, not a FASTA file, no header
#' @param motif_file assumed to be in HOMER format
#' @param outfile where score matrix is written
#' @param HOMER_PATH path to HOMER bin directory
create_motif_score_matrix <- function(sequence_file,
                                      motif_file,
                                      outfile,
                                      ncore,
                                      HOMER_PATH)
{
  seqs <- data.table::fread(sequence_file, header=F) %>% as.matrix %>% as.character()
  dna <- DNAStringSet(seqs) %>% setNames(paste("seq", 1:length(seqs), sep=""))
  control <- dna[1] %>% setNames("control_seq")
  
  seq_fasta <- paste("temp_sequences", round(1E10*runif(1)), ".fasta", sep="")
  control_fasta <- paste("temp_control", round(1E10*runif(1)), ".fasta", sep="")
  writeXStringSet(dna, seq_fasta)
  writeXStringSet(control, control_fasta)
  
  temp_output_folder <- paste("temp_output", round(1E10*runif(1)), "/", sep="")
  dir.create(temp_output_folder)
  
  motifs <- read_homer_motifs.BM(motif_file)
  names(motifs) <- sapply(names(motifs), function(m) gsub("::", "_", m))
  ptm <- proc.time()
  
  mclapply(names(motifs), function(mn) {
    cat("calling homer with motif", mn, "\n")
    tb <- find_motifs.homer(seq_fasta, 
                            control_fasta, 
                            motifs[mn], 
                            just_plus_strand=F, 
                            HOMER_PATH)
    cat("finding max scores for motif", mn, "\n")
    tb2 <- HOMER2MaxScore.homer(tb, names(dna))
    cat("writing max scores for motif", mn, "\n")
    temp_outfile <- paste(temp_output_folder, mn, "_matrix_scores.csv", sep="")
    write.csv(tb2, temp_outfile, row.names=F)
  }, mc.cores=ncore)
  
  all_output_files <- dir(temp_output_folder)
  output_motifs <- sapply(strsplit(all_output_files, split="_matrix"), "[", 1)
  m <- sapply(all_output_files, function(f) {
    read.csv(paste(temp_output_folder, f, sep=""), 
             header=T, stringsAsFactors = F) %>% dplyr::pull(score)
  })
  colnames(m) <- output_motifs
  
  write.csv(m, outfile, row.names = F)
  file.remove(seq_fasta)
  file.remove(control_fasta)
  Sys.setenv(tf=temp_output_folder)
  system("rm -r $tf")
  
  atm <- proc.time()
  print(atm - ptm)
}



