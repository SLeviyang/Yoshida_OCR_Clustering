library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)

#' Find all ISRE loci on the genome
create_all_ISRE_loci_file <- function(species="mouse")
{
  all_motifs <- get_ISRE_search_motifs()

  cat("loading the", species, "genome\n")
  g <- BSgenome.Mmusculus.UCSC.mm10
  #g <- GenomeR::get.genome(species)
  include_chr <- intersect(names(g), c("chrX", paste("chr", 1:100, sep="")))
  target_seqs <- BSgenome::getSeq(g, names=include_chr)

  all_tb <- plyr::adply(all_motifs, 1, function(cmotif) {
    cat("searching for for base motif", cmotif, "\n")

    tb <- motif2track(motif=cmotif, target_seqs=target_seqs)
    # put in BED6 format
    tb <- dplyr::mutate(tb, peak_name=paste("peak", 1:nrow(tb), "_", cmotif, sep=""), score=".") %>%
      dplyr::select(chr, chrStart, chrEnd, peak_name, score, strand, motif, revC)

    return (tb)
  }, .id=NULL)

  data.table::fwrite(all_tb, "../../all_ISRE.bed", sep="\t", col.names = F)
  system("bedtools sort -i ../../all_ISRE.bed > ../../all_ISRE_sorted.bed")
}

get_ISRE_search_motifs <- function()
{
  acgt <- c("A", "C", "G", "T")
  wild_indices <- c(1:4, 7:10)
  base_motif <- "TTTCNNTTTC"
  base_motif_split <- strsplit(base_motif, split="")[[1]]
  
  all_motifs <- sapply(wild_indices, function(i) {
    sapply(acgt, function(x) {
      s <- base_motif_split
      s[i] <- x
      return (paste(s, collapse=""))
    })
  }) %>% as.character()
  all_motifs <- c(base_motif, all_motifs) %>% unique
  
  return (all_motifs)
}



#' Find motifs within sequences
#'
#' @param motif a string containing the motif to be searched for
#' @param target_seqs a DNAStringSet giving the targets over which
#' motif is searched for.  Typically a collection of chromosomes.
#'
#' @details the reverse complement of the motif is also searched for.
#'
#' @returns a data.frame with fields chr, chrStart, chrEnd, strand, motif, revC.
#' revC is true if match motif is the reverse complement of the passed motif.  The
#' motif field contains the match motif, which may be the reverse complement.
motif2track <- function(motif, target_seqs)
{
  exploded_motifs <- lapply(strsplit(motif, split="")[[1]], function(n)
    if (n=="N")
      return (c("A", "C", "G", "T"))
    else
      return (n)
  ) %>%
    expand.grid(stringsAsFactors = F) %>%
    as.matrix %>%
    apply(1, paste, collapse="") %>%
    Biostrings::DNAStringSet()
  revC_exploded_motifs <- Biostrings::reverseComplement(exploded_motifs)

  # I look for motifs only on the positive strand!  I'll store whether I found
  # the motif passed or its reverse complement
  motif_tb <- rbind(data.frame(motif=as.character(exploded_motifs), revComplement=F, stringsAsFactors = F),
                    data.frame(motif=as.character(revC_exploded_motifs), revComplement=T, stringsAsFactors = F))

  out <- plyr::adply(1:length(target_seqs), 1, function(i) {
    cat("pattern matching for", names(target_seqs)[i], "\n")
    x <- plyr::adply(motif_tb, 1, function(tb) {
      m <- tb$motif; rc <- tb$revComplement
      z <- Biostrings::vmatchPattern(pattern=m, subject=target_seqs[i], fixed=T)
      as.data.frame(z[[1]]) %>% dplyr::mutate(chr=names(target_seqs)[i],
                                              motif=m,
                                              strand=1,
                                              revC=rc)
    })
    return (x)
  }, .id=NULL) %>%
    dplyr::rename(chrStart=start, chrEnd=end) %>%
    dplyr::select(chr, chrStart, chrEnd, strand, motif, revC)

  return (out)
}


################################
create_all_ISRE_loci_file()

