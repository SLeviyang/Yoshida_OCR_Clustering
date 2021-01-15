load_JASPAR_motifs <- function()
{
  opt <- list(tax_group="vertebrates", matrixtype="PWM")
  p <- getMatrixSet(JASPAR2018, opt)
  return (p)
}


# read in chromVar mouse v2 cisDB motifs used by Yoshida et al as TF motifs
# downloaded from
# https://github.com/buenrostrolab/chromVARmotifs/blob/master/data/mouse_pwms_v2.rda
# Note:  Yoshida Table S5 contains a subset of 430 of the version 2 motifs, but these
# 430 do not include SPI-1.   
load_chromVARmotifs_motifs <- function(version=2)
{
  if (version==1) {
    load("../../data_input/mouse_pwms_v1.rda")
    motifs <- mouse_pwms_v1
  } else {
    load("../../data_input/mouse_pwms_v2.rda")
    motifs <- mouse_pwms_v2
  }
  motif_names <- sapply(motifs, function(m) name(m)) %>% setNames(NULL)
  names(motifs) <- toupper(motif_names)
  
  return (motifs)
}

create_motif_score_matrices <- function(max_cluster_number=14)
{
  motifCluster_dir <- paste(get_DATA_DIR(), "TF_motifs/", sep="")
  if (!dir.exists(motifCluster_dir))
    dir.create(motifCluster_dir)

  motifs <- load_chromVARmotifs_motifs(version=2)
  motif_names <- sapply(motifs, name) %>% toupper
  if (length(motif_names) != length(unique(motif_names)))
    stop("motif names are not unique!")
  pc <- construct.peak_clusters()
  
  mclapply(0:max_cluster_number, function(cluster_number) {
    motif_score_file <- paste(motifCluster_dir, "motif_scores_cluster_", 
                              cluster_number, ".csv", sep="")
    rmotif_score_file <- paste(motifCluster_dir, "motif_scores_cluster_", 
                              cluster_number, "_permuted.csv", sep="")
    print(motif_score_file)
    print(rmotif_score_file)
    
    seqs <- form_cluster_sequences.peak_clusters(pc, cluster_number)
    
    nc <- nchar(seqs[1])
    rseqs <- sapply(strsplit(seqs, split=""), function(s) 
                    s[sample.int(nc)] %>% paste(collapse=""))
    
    cc <- motifmatchr::matchMotifs(motifs, seqs, out="scores", p.cutoff=1) 
    rcc <- motifmatchr::matchMotifs(motifs, rseqs, out="scores", p.cutoff=1) 
    
    m <-  SummarizedExperiment::assay(cc) %>% as.matrix
    colnames(m) <- motif_names
    rm <-  SummarizedExperiment::assay(rcc) %>% as.matrix
    colnames(rm) <- motif_names
   
    write.csv(m, motif_score_file, row.names = F)
    write.csv(rm, rmotif_score_file, row.names = F)
  }, mc.cores = 5)
  
}

#################################################
construct.TF <- function(max_cluster_number=14)
{
  motifCluster_dir <- paste(get_DATA_DIR(), "TF_motifs/", sep="")
  all_clusters <- 0:max_cluster_number

  m_list <- vector(mode="list", length=length(all_clusters))
  rm_list <- vector(mode="list", length=length(all_clusters))
  cluster_list <- vector(mode="list", length=length(all_clusters))

  for (i in 1:length(all_clusters)) {
    cluster_number <- all_clusters[i]
    cat("loading motif scores for cluster", cluster_number, "\n")
    file <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, ".csv", sep="")
    rfile <- paste(motifCluster_dir, "motif_scores_cluster_", cluster_number, "_permuted.csv", sep="")

    m_list[[i]] <- data.table::fread(file, sep=",", header=T) %>% as.matrix
    rm_list[[i]] <- data.table::fread(rfile, sep=",", header=T) %>% as.matrix
    cluster_list[[i]] <- rep(cluster_number, nrow(m_list[[i]]))
  }

  m <- do.call(rbind, m_list)
  rm <- do.call(rbind, rm_list)
  cluster <- do.call(c, cluster_list)

  return (list(scores=m_list,
               null_scores=rm_list))
}


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
fit_within.TF <- function(mc, index1, index2, ncv=5, verbose=T)
{
  scores <- mc$scores[[index1+1]]
  null_scores <- mc$scores[[index2+1]]
  
  r <- c(rep(1, nrow(scores)), rep(0, nrow(null_scores)))
  X <- rbind(scores, null_scores) %>% scale
  browser()

  if (verbose)
    cat("number samples=", nrow(X), "nfeatures=", ncol(X),"\n")


  foldid <- rep(0, length(r))
  foldid[r==1] <- rep(1:ncv, ceiling(sum(r==1)/ncv))[1:sum(r==1)]
  foldid[r==0] <- rep(1:ncv, ceiling(sum(r==0)/ncv))[1:sum(r==0)]

  g_cv <- cv.glmnet(X, r, family="binomial", foldid=foldid)
  g <- g_cv$glmnet.fit

  optimal_ind <- which(g_cv$lambda == g_cv$lambda.1se)
  coefs <- (g$beta %>% as.matrix)[,optimal_ind]
  a0 <- g$a0[optimal_ind]
  browser()

  names(a0) <- "a0"

  return (list(a0=a0, coefs=coefs))
}