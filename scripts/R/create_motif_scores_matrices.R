load_JASPAR_motifs <- function()
{
  meme_file <- paste(get_INPUT_DATA_DIR(),
                     "JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt",
                     sep="")
  homer_file <- paste("temp", round(1E10*runif(1)), ".txt", sep="")
  meme2homer.BM(meme_file, homer_file, min_nsites=0, set_score=-1E10)
  
  motifs <- read_homer_motifs.BM(homer_file)
  file.remove(homer_file)
  
  full_motif_names <- get_motif_names.BM(motifs) %>% toupper
  motif_names <- sapply(strsplit(full_motif_names, split="_"), "[", 2)
 
  if (length(motif_names) != length(unique(motif_names)))
    stop("motif names are not unique!")
  names(motifs) <- motif_names
  
  return (motifs)
}

# load_JASPAR_motifs <- function()
# {
#   opt <- list(tax_group="vertebrates", matrixtype="PWM")
#   p <- getMatrixSet(JASPAR2018, opt)
#   return (p)
# }
# 
# 
# # read in chromVar mouse v2 cisDB motifs used by Yoshida et al as TF motifs
# # downloaded from
# # https://github.com/buenrostrolab/chromVARmotifs/blob/master/data/mouse_pwms_v2.rda
# # Note:  Yoshida Table S5 contains a subset of 430 of the version 2 motifs, but these
# # 430 do not include SPI-1.   
# load_chromVARmotifs_motifs <- function(version=2)
# {
#   if (version==1) {
#     load("../../data_input/mouse_pwms_v1.rda")
#     motifs <- mouse_pwms_v1
#   } else {
#     load("../../data_input/mouse_pwms_v2.rda")
#     motifs <- mouse_pwms_v2
#   }
#   motif_names <- sapply(motifs, function(m) name(m)) %>% setNames(NULL)
#   names(motifs) <- toupper(motif_names)
#   
#   return (motifs)
# }

create_all_motif_score_matrices <- function(max_cluster_number=20)
{
  motifs <- load_JASPAR_motifs()
  create_motif_score_matrices(max_cluster_number, motifs, distal=T)
  create_motif_score_matrices(max_cluster_number, motifs, distal=F)
}

create_motif_score_matrices <- function(max_cluster_number,
                                        motifs,
                                        distal=TRUE)
{
  motifCluster_dir <- paste(get_DATA_DIR(), "TF_motifs/", sep="")
  if (!dir.exists(motifCluster_dir))
    dir.create(motifCluster_dir)
  if (distal)
    raw_motifCluster_dir <- paste(motifCluster_dir, "raw_matrices_distal/", sep="")
  else
    raw_motifCluster_dir <- paste(motifCluster_dir, "raw_matrices_proximal/", sep="")
  if (!dir.exists(raw_motifCluster_dir))
    dir.create(raw_motifCluster_dir)

  pc <- construct.peak_clusters()
  
  mclapply(0:max_cluster_number, function(cluster_number) {
    motif_score_file <- paste(raw_motifCluster_dir, "motif_scores_cluster_", 
                              cluster_number, ".csv", sep="")
    rmotif_score_file <- paste(raw_motifCluster_dir, "motif_scores_cluster_", 
                              cluster_number, "_permuted.csv", sep="")
    sequence_score_file <- paste(raw_motifCluster_dir, "motif_scores_cluster_", 
                              cluster_number, "_sequences.csv", sep="")
    print(motif_score_file)
    print(rmotif_score_file)
    
    seqs <- form_cluster_sequences.peak_clusters(pc, cluster_number)
    g <- form_cluster_genomics.peak_clusters(pc, cluster_number)
    
    
    if (distal)
      active_loci <- g$dTSS > 3000
    else
      active_loci <- g$dTSS < 500
    
  
    seqs <- seqs[active_loci]
    g <- g[active_loci,]
    
    nc <- nchar(seqs[1])
    rseqs <- sapply(strsplit(seqs, split=""), function(s) 
                    s[sample.int(nc)] %>% paste(collapse=""))
    
    if (!file.exists(motif_score_file)) {
      m <- HOMERscoreMatrix(seqs, motifs, just_plus_strand = F,
                            message=paste("cluster", cluster_number, sep=" "))
      write.csv(m, motif_score_file, row.names = F)
    }
    if (!file.exists(rmotif_score_file)) {
      rm <-  HOMERscoreMatrix(rseqs, motifs, just_plus_strand = F,
                          message=paste("null cluster", cluster_number, sep=" "))
      write.csv(rm, rmotif_score_file, row.names = F)
    }
    if (!file.exists(sequence_score_file)) {
      write.csv(data.frame(sequence=seqs, stringsAsFactors = F), 
                sequence_score_file, row.names = F)
    }
  }, mc.cores = 3)
  
}

#################################################
load_motif_score_matrix <- function(cluster_number, permuted=F)
{
  motifCluster_dir <- paste(get_DATA_DIR(), "TF_motifs/", sep="")
  cat("loading motif scores for cluster", cluster_number, "\n")
  
  if (!permuted)
    file <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, ".csv", sep="")
  else 
    file <- rfile <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, 
                           "_permuted.csv", sep="")
  
  return (data.table::fread(file, sep=",", header=T) %>% as.matrix)
}
# construct.TF <- function(max_cluster_number=14)
# {
#   motifCluster_dir <- paste(get_DATA_DIR(), "TF_motifs/", sep="")
#   all_clusters <- 0:max_cluster_number
# 
#   m_list <- vector(mode="list", length=length(all_clusters))
#   rm_list <- vector(mode="list", length=length(all_clusters))
#   cluster_list <- vector(mode="list", length=length(all_clusters))
# 
#   for (i in 1:length(all_clusters)) {
#     cluster_number <- all_clusters[i]
#     cat("loading motif scores for cluster", cluster_number, "\n")
#     file <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, ".csv", sep="")
#     rfile <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, "_permuted.csv", sep="")
# 
#     m_list[[i]] <- data.table::fread(file, sep=",", header=T) %>% as.matrix
#     rm_list[[i]] <- data.table::fread(rfile, sep=",", header=T) %>% as.matrix
#     cluster_list[[i]] <- rep(cluster_number, nrow(m_list[[i]]))
#   }
# 
#   m <- do.call(rbind, m_list)
#   rm <- do.call(rbind, rm_list)
#   cluster <- do.call(c, cluster_list)
# 
#   return (list(scores=m_list,
#                null_scores=rm_list))
# }


# ###########################################################################
# form_cluster_motif_matrix.motifCluster <- function(mc, index)
# {
#   m <- mc$scores[mc$cluster == index,]
#   return (m)
# }
# 
# show_motif_scores.motifCluster <- function(mc, motif, 
#                                            nrow, ncol,
#                                            clusters=NULL)
# {
#   d <- rbind(data.frame(score=mc$scores[,motif], cluster=mc$cluster, type="data"),
#              data.frame(score=mc$null_scores[,motif], cluster=mc$cluster, type="null"))
#   
#   if (!is.null(clusters))
#     d <- dplyr::filter(is.element(cluster, clusters))
#   
#   g <- ggplot(data=d) + stat_ecdf(aes(x=score, color=type)) + 
#     facet_wrap("cluster", nrow=nrow, ncol=ncol)
#   
#   return (g)
# }
# 
# fit_within.TF <- function(mc, index1, index2, ncv=5, verbose=T)
# {
#   scores <- mc$scores[[index1+1]]
#   null_scores <- mc$scores[[index2+1]]
#   
#   r <- c(rep(1, nrow(scores)), rep(0, nrow(null_scores)))
#   X <- rbind(scores, null_scores) %>% scale
#   browser()
# 
#   if (verbose)
#     cat("number samples=", nrow(X), "nfeatures=", ncol(X),"\n")
# 
# 
#   foldid <- rep(0, length(r))
#   foldid[r==1] <- rep(1:ncv, ceiling(sum(r==1)/ncv))[1:sum(r==1)]
#   foldid[r==0] <- rep(1:ncv, ceiling(sum(r==0)/ncv))[1:sum(r==0)]
# 
#   g_cv <- cv.glmnet(X, r, family="binomial", foldid=foldid)
#   g <- g_cv$glmnet.fit
# 
#   optimal_ind <- which(g_cv$lambda == g_cv$lambda.1se)
#   coefs <- (g$beta %>% as.matrix)[,optimal_ind]
#   a0 <- g$a0[optimal_ind]
#   browser()
# 
#   names(a0) <- "a0"
# 
#   return (list(a0=a0, coefs=coefs))
# }