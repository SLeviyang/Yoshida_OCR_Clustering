# An R version of the peak_clusters class defined in peak_clusters.py

construct.peak_clusters <- function(FDR_string="3.0")
{
  # DATA_DIR defined in configuration.R
  cluster_file <- paste(get_DATA_DIR(), "peak_clusters/clusters_",
                        FDR_string, ".csv", sep="")
  sequence_file <- paste(get_DATA_DIR(), "master_peaks/master_peak_sequences.txt", sep="")
  matrix_file <- paste(get_DATA_DIR(), "master_peaks/master_peak_matrix.csv", sep="")
  genomics_file =  paste(get_DATA_DIR(), "master_peaks/master_peak_genomic_info.csv", sep="")
  
  pc <- list(clusters=read.csv(cluster_file, header=T, stringsAsFactors = F),
             sequences=read.csv(sequence_file, header=F, stringsAsFactors = F),
             matrix=data.table::fread(matrix_file, sep=",", header=T),
             genomics=data.table::fread(genomics_file, sep=",", header=T, stringsAsFactors = F))
  
  return (pc)
}



form_cluster_sequences.peak_clusters <- function(pc, index, window_size=100)
{
  cl <- pc$clusters
  rows <- dplyr::filter(cl, cluster == index) %>%
          dplyr::pull(row)
  seqs <- pc$sequences[rows+1,]
  if (length(nchar(seqs) %>% unique) != 1)
    stop("sequences should all be the same length!")
  nc <- nchar(seqs[1])
  seqs_m <- matrix(strsplit(seqs, split="") %>% unlist, byrow=T, ncol=nchar(seqs[1]))
  
  mid_point <- round(nc/2)
  start_col = mid_point - window_size
  end_col = mid_point + window_size
  if (start_col > nc | end_col < 1)
    stop("window is too wide!")
  
  seqs_windowed <- apply(seqs_m[,start_col:end_col], 1, paste, collapse="")
  
  return (seqs_windowed)
}

form_cluster_matrix.peak_clusters <- function(pc, index)
{
  cl <- pc$clusters
  rows <- dplyr::filter(cl, cluster == index) %>%
    dplyr::pull(row)
  m <- pc$matrix[rows+1,]
  
  return (m)
}

form_cluster_genomics.peak_clusters <- function(pc, index)
{
  cl = pc$clusters
  rows <- dplyr::filter(cl, cluster == index) %>%
    dplyr::pull(row)
  g = pc$genomics[rows+1,]
  
  return (g)
}

