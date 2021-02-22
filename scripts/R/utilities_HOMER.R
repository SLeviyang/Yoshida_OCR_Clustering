# functions to search/define biniding motifs

#' Use HOMER utility to find binding motifs in sequences
#'
#' @param seq_file fasta file
#' @param control_file fasta
#' @param motif_list a list of PWM
#' http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html for format)
#' @param just_plus_strand should the reverse complement be searched?
#' @param HOMER_BIN_dir path to HOMER's bin/
#' @score_matrix should a score matrix be returned, or the raw HOMER output?
#'
#' @details if both seqs and fasta file are not NULL, seqs is used
#'
#' @return a data.frame containing HOMER output
find_motifs.homer <- function(seq_file,
                             control_file,
                             motif_list,
                             just_plus_strand,
                             HOMER_dir,
                             use_cumulative_binomial=T)
{
  if (length(motif_list) != 1)
    return ("only one motif is allowed, for more just call function repeatedly")
  
  motif_file <- paste("temp_homer_motif", round(1E6*runif(1)), ".txt", sep="")
  write_homer_motifs.BM(motif_list, motif_file)
  nmotifs <- length(motif_list)
  
  cat("processing", nmotifs, "motifs\n")

  out_dir <- paste("homer_workdir_temp", round(1E8*runif(1)), "/", sep="")
  dir.create(out_dir)
  outfile <- paste(out_dir, "HOMER_analysis.txt", sep="")

  spath <- Sys.getenv("PATH")
  # check if homer is already there
  if (!grepl(HOMER_dir, spath))
    spath <- paste(spath, HOMER_dir, sep=":")

  Sys.setenv(f1=seq_file, f2=control_file, mf=motif_file,
             od=out_dir, PATH=spath,
             of=outfile)
  if (just_plus_strand)
    Sys.setenv(ps="-norevopp")
  else
    Sys.setenv(ps="")
  if (use_cumulative_binomial)
    Sys.setenv(cb="-b")
  else
    Sys.setenv(cb="")
  system("findMotifs.pl $f1 fasta $od -fastaBg $f2 -find $mf $cb  $ps > $of")

  # in return file, offset is relative to center nucleatide

  #tb <- read.delim(outfile, header=T, stringsAsFactors = F)
  print("reading homer output file into R...")
  tb <- data.table::fread(outfile, header=T, stringsAsFactors = F, sep="\t") %>%
        setNames(c("FASTA.ID", "Offset", "Sequence", "Motif.Name", 
                   "strand", "MotifScore"))

  for (f in dir(out_dir, full.names = T))
    file.remove(f)
  system("rm -r $od")
  file.remove(motif_file)

  return  (tb)
}


#' Get the max score of a motif for each sequence using homer
#' 
#' @param tb a data.frame returned by find_motifs.homer
#' @param seq_names names of sequences in fasta files sent to find_motifs.homer
HOMER2MaxScore.homer <- function(homer_tb, seq_names)
{
  score_tb <- plyr::ddply(homer_tb, .(FASTA.ID), function(ctb) {
    data.frame(score=max(ctb$MotifScore), stringsAsFactors = F)
  }) %>% dplyr::rename(id=FASTA.ID)
  missing_seqs <- setdiff(seq_names, score_tb$id)
  if (length(missing_seqs) > 0) 
    score_tb <- rbind(score_tb,
                      data.frame(id=missing_seqs, score=NA))
  
  score_tb <- score_tb[match(seq_names, score_tb$id),]
  return (score_tb)
}


# A BM object is a motif list.
# Functions to read and write MEME and HOMER motif files
# A motif list object is a list of motifs.  Each
# motif is a list with fields motif_name, sequence, score, nsites, m

#' Get the motif names in a motif list
get_motif_names.BM <- function(ml)
{
  return (sapply(ml, "[[", "motif_name"))
}

get_sequences.BM <- function(ml)
{
  return (sapply(ml, "[[", "sequence"))
}

# methods for HOMER motif format

#' Convert homer motif file to homer object
#'
#' @param homer_file a path to a homer motif file
#'
#' @details a homer_motifs object is a list of
#' motif objects.  A motif object is a list with
#' fields m (matrix), seq (conensus_seq), motif_name.
#'
#' Assumes that homer motif files are tab delimited
#'
#' @return a motif list
read_homer_motifs.BM <- function(homer_file)
{
  f <- file(homer_file, open = "r")
  rl <- readLines(f)
  # remove all empty lines
  rl <- rl[rl != ""]

  start <- grep(">", rl)
  end <- c(start[-1]-1, length(rl))
  n_motifs <- length(start)

  hm <- NULL
  for (i in 1:n_motifs) {
    s <- strsplit(rl[start[i]], split="\t")[[1]]
    sequence <- substr(s[1], 2, nchar(s[1]))
    motif_name <- s[2]
    score <- s[3] %>% as.numeric

    m <- NULL
    for (j in (start[i]+1):end[i]) {
      s <- strsplit(rl[j], split="\t")[[1]]
      m <- rbind(m, s %>% as.numeric)
    }

    # check
    if (nrow(m) != nchar(sequence))
      stop("error reading homer motif file")

    hm <- c(hm,
            list(list(sequence=sequence,
                      motif_name=motif_name,
                      score=score,
                      nsites=NA,
                      m=m)))
  }
  cat("number of motifs read", length(hm), "\n")
  close(f)

  names(hm) <- sapply(hm, "[[", "motif_name")

  return (hm)
}

write_homer_motifs.BM <- function(motif_list, out_homer_file)
{
  # now convert to HOMER format
  f <- file(out_homer_file, open = "w")
  for (info in motif_list) {
    s <- paste(">", info$sequence, "\t", info$motif_name, "\t", info$score, sep="")
    writeLines(s, f)
    apply(info$m, 1, function(r) writeLines(paste(r, collapse="\t"), f))
    writeLines("", f)
  }

  cat("number of motifs written", length(motif_list), "\n")
  close(f)
}

# methods for MEME motif format

#' Parses the meme.txt output file for PWM
#'
#' @param f path to meme.txt output file
#'
#' @return a motif object (A motif object is a list with
#' fields m (matrix), seq (conensus_seq), motif_name.) and
#' a data.frame containing fields gene,
parse_output.meme <- function(meme_motif_file, set_score=NA)
{
  f <- file(meme_motif_file, open="r")
  rl <- readLines(f)

  motif_lines <- grep("letter-probability matrix", rl)

  get_motif_specification <- function(motif_info_parsed, char) {
    ind <- which(motif_info_parsed==char)
    if (length(ind)==0)
      return (NA)
    else
      return (motif_info_parsed[ind+1] %>% as.numeric)
  }

  ACGT <- c("A", "C", "G", "T")
  motif_list <- plyr::alply(motif_lines, 1, function(start) {

    # example of motif name line
    #  MOTIF cr
    motif_name <- strsplit(rl[start-2], split="[\t ]")[[1]][3]

    # example of motif info line
    # letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
    motif_info_parsed <- strsplit(rl[start], split="[ =]")[[1]]
    motif_info_parsed <- motif_info_parsed[motif_info_parsed != ""]

    w <- get_motif_specification(motif_info_parsed, "w")
    nsites <- get_motif_specification(motif_info_parsed, "nsites")

    m <- matrix(0, nrow=w, ncol=4)
    seq <- vector("character", length=w)
    for (i in 1:w) {
      vals <- strsplit(rl[start+i], split=" ")[[1]]
      vals <- vals[vals != ""] %>% as.numeric
      m[i,] <- vals
      seq[i] <- ACGT[which.max(m[i,])]
    }


    info <- list(motif_name=motif_name,
                 sequence=paste(seq, collapse=""),
                 m=m,
                 score=set_score,
                 nsites=nsites)
    return (info)
  })

  names(motif_list) <- sapply(motif_list, function(m) m$motif_name)

  # now get matches
  nsites <- motif_list[[1]]$nsites
  match_line <- grep("MEME-1 in BLOCKS format", rl)
  cols <- c("sequence_name",  "start",  "motif", "unknown")
  m <- matrix("", nrow=nsites, ncol=length(cols))
  for (i in 1:nsites) {
    s <- strsplit(rl[match_line+2+i], split="[ ()]")[[1]]
    s <- s[s != ""]
    m[i,] <- s
  }
  tb <- base::as.data.frame(m, stringsAsFactors=F) %>%
          setNames(cols) %>%
          dplyr::select(sequence_name, start, motif)
  tb$start <- as.numeric(tb$start)

  close(f)
  return (list(motif=motif_list[[1]], tb=tb))
}

# Converts meme motifs to a motif list
read_meme_motifs.BM <- function(meme_motif_file, set_score=NA)
{
  f <- file(meme_motif_file, open="r")
  rl <- readLines(f)

  motif_lines <- grep("MOTIF", rl)

  get_motif_specification <- function(motif_info_parsed, char) {
    ind <- which(motif_info_parsed==char)
    if (length(ind)==0)
      return (NA)
    else
      return (motif_info_parsed[ind+1] %>% as.numeric)
  }

  ACGT <- c("A", "C", "G", "T")
  motif_list <- plyr::alply(motif_lines, 1, function(start) {

    # example of motif name line
    #  MOTIF cr
    motif_name <- strsplit(rl[start], split=" ")[[1]][-1] %>%
      paste(collapse="_")

    # example of motif info line
    # letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
    motif_info_parsed <- strsplit(rl[start+1], split="[ =]")[[1]]
    motif_info_parsed <- motif_info_parsed[motif_info_parsed != ""]

    w <- get_motif_specification(motif_info_parsed, "w")
    nsites <- get_motif_specification(motif_info_parsed, "nsites")

    m <- matrix(0, nrow=w, ncol=4)
    seq <- vector("character", length=w)
    for (i in 1:w) {
      vals <- strsplit(rl[start+1+i], split=" ")[[1]]
      vals <- vals[vals != ""] %>% as.numeric
      m[i,] <- vals
      seq[i] <- ACGT[which.max(m[i,])]
    }


    info <- list(motif_name=motif_name,
                 sequence=paste(seq, collapse=""),
                 m=m,
                 score=set_score,
                 nsites=nsites)
    return (info)
  })

  close(f)
  return (motif_list)
}

# motif file format conversion methods 
#' Convert a MEME motif file to a HOMER motif file
#'
#' @param min_nsites nsites gives the number of ChIPseq sites on
#' which the data is based and is given in MEME files, but not homer files
#' @param set_score score is used in HOMER as the cutoff log likelihood for 
#' scoring a subsequence of a sequence.  A -Inf will lead to every position being scored.
meme2homer.BM <- function(meme_file, out_homer_file,
                          min_nsites=0,
                          set_score=-1000000)
{
  motif_list <- read_meme_motifs.BM(meme_file, set_score=set_score)
  motif_list <- motif_list[sapply(motif_list, function(mot) mot$nsites >= min_nsites)]

  write_homer_motifs.BM(motif_list, out_homer_file)
}

