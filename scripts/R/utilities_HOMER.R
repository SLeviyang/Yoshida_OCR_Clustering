# functions to search/define biniding motifs

#' Use HOMER utility to find binding motifs in sequences
#'
#' @param seqs A DNAStringSet
#' @param fasta_file path to a fasta file
#' @param motif_file path to a HOMER motif file (see
#' http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html for format)
#' @param just_plus_strand should the reverse complement be searched?
#' @param HOMER_BIN_dir path to HOMER's bin/
#'
#' @details if both seqs and fasta file are not NULL, seqs is used
#'
#' @return a data.frame containing HOMER output
find_motifs.homer <- function(seqs,
                             control_seqs,
                             motif_file,
                             just_plus_strand=F,
                             HOMER_dir=HOMER_PATH)
{

  out_dir <- paste("homer_workdir_temp", round(1E8*runif(1)), "/", sep="")
  dir.create(out_dir)
  outfile <- paste(out_dir, "HOMER_analysis.txt", sep="")
  f1 <- paste("seqs_HOMMER", round(1E8*runif(1)), ".fasta", sep="")
  f2 <- paste("control_HOMMER", round(1E8*runif(1)), ".fasta", sep="")

  writeXStringSet(seqs, f1)
  writeXStringSet(control_seqs, f2)

  spath <- Sys.getenv("PATH")
  # check if homer is already there
  if (!grepl(HOMER_dir, spath))
    spath <- paste(spath, HOMER_dir, sep=":")

  Sys.setenv(f1=f1, f2=f2, mf=motif_file,
             od=out_dir, PATH=spath,
             of=outfile)
  if (just_plus_strand)
    Sys.setenv(ps="-norevopp")
  else
    Sys.setenv(ps="")
  system("findMotifs.pl $f1 fasta $od -fastaBg $f2 -find $mf $ps > $of")

  # in return file, offset is relative to center nucleatide

  #tb <- read.delim(outfile, header=T, stringsAsFactors = F)
  print("reading homer output file into R...")
  tb <- data.table::fread(outfile, header=T, stringsAsFactors = F, sep="\t") %>%
        setNames(c("FASTA.ID", "Offset", "Sequence", "Motif.Name", 
                   "strand", "MotifScore"))

  tb$strand <- ifelse(tb$strand=="+", 1, -1)
  system("rm $f1")
  system("rm $f2")

  for (f in dir(out_dir, full.names = T))
    file.remove(f)
  system("rm -r $od")

  return  (tb)
}

#' Run homer findMotifs.pl on DNAStringSet
#'
#' @param gene_sequences DNAStringSet
#' @param control_sequences DNAStringSet with no effect
#' @param motif_list a list of PWM
#' @param just_plus_strand should just the + strand be used
#'
#' @details HOMER doesn't seem to use the control sequences for scoring
#'
#' @return a data.frame containing HOMER output
sequences2homer <- function(gene_sequences,
                            control_sequences,
                            motif_list,
                            just_plus_strand)
{
  temp_motif_file <- paste("temp_homer_motif", round(1E6*runif(1)), ".txt", sep="")
  write_homer_motifs.BM(motif_list, temp_motif_file)
  nmotifs <- length(motif_list)

  cat("processing", nmotifs, "motifs\n")
  all_seqs <- c(gene_sequences, control_sequences)
  homer <- find_motifs.homer(gene_sequences,
                    control_sequences,
                    temp_motif_file,
                    just_plus_strand) %>%
    dplyr::rename(id=FASTA.ID)

  Sys.setenv(tf=temp_motif_file)
  system("rm $tf")

  return (homer)
}

#' Create a score matrix (seqs by motifs)
#' 
#' @param seqs a character.vector
#' @param motifs a names list of BM objects, see below
HOMERscoreMatrix <- function(seqs, motifs, just_plus_strand=F,
                             message=NA)
{
  motif_names <- names(motifs)
  nm <- length(motifs)
  dna_seqs <- DNAStringSet(seqs) %>% 
              setNames(paste("seq", 1:length(seqs), sep=""))
  
  m <- sapply(1:length(motifs), function(i) {
    if (!is.na(message))
      print (message)
    cat("processing motif", i, "of", nm, ":", motif_names[i], "\n")
    tb <- sequences2homer(dna_seqs, dna_seqs[1], motifs[i], just_plus_strand)
    # when a probability in the PWM is 0, HOMER might not return a score
  
    missing_seqs <- setdiff(names(dna_seqs), tb$id)
    if (length(missing_seqs) > 0) {
      tb_missing <- plyr::adply(missing_seqs, 1, function(sn) {
        data.frame(id=sn, Offset=NA, Sequence=NA,
                   Motif.Name=motif_names[i], strand=NA, MotifScore=NA)
      }, .id=NULL)
      tb <- rbind(tb, tb_missing)
    }
    plyr::daply(tb, .(id), function(tbm) {
      if (all(is.na(tbm$MotifScore)))
        return (NA)
      else
        return (max(tbm$MotifScore, na.rm=T))
    })
  })
  if (length(dna_seqs)==1) {
    m <- matrix(m, nrow=1)
  } else if (length(motif_names)==1) {
    m <- matrix(m, ncol=1)
  }
  colnames(m) <- motif_names
  
  return (m)
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

