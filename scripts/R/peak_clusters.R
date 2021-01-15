# An R version of the peak_clusters class defined in peak_clusters.py

construct.peak_clusters <- function()
{
  # DATA_DIR defined in configuration.R
  cluster_file <- paste(get_DATA_DIR(), "peak_clusters/clusters.csv", sep="")
  sequence_file <- paste(get_DATA_DIR(), "master_peaks/master_peak_sequences.txt", sep="")
  matrix_file <- paste(get_DATA_DIR(), "master_peaks/master_peak_matrix.csv", sep="")
  
  pc <- list(clusters=read.csv(cluster_file, header=T, stringsAsFactors = F),
             sequences=read.csv(sequence_file, header=F, stringsAsFactors = F),
             matrix=data.table::fread(matrix_file, sep=",", header=T))
  
  return (pc)
}



form_cluster_sequences.peak_clusters <- function(pc, index)
{
  cl <- pc$clusters
  rows <- dplyr::filter(cl, cluster == index) %>%
          dplyr::pull(row)
  seqs <- pc$sequences[rows+1,]
  
  return (seqs)
}

form_cluster_matrix.peak_clusters <- function(pc, index)
{
  cl <- pc$clusters
  rows <- dplyr::filter(cl, cluster == index) %>%
    dplyr::pull(row)
  m <- pc$matrix[rows+1,]
  
  return (m)
}

