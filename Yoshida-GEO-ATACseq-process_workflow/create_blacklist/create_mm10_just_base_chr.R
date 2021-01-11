library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)

# first create fasta files
g <- BSgenome.Mmusculus.UCSC.mm10
base_names <- c(paste("chr", 1:19, sep=""), "chrX")
base_fa <- "mm10_base_chromosomes.fa"
if (!file.exists(base_fa)) {
  seq <- getSeq(g, base_names)
  writeXStringSet(seq, base_fa)
}

if (!file.exists("chrM.fa")) {
  chrM <- getSeq(g, "chrM")
  writeXStringSet(DNAStringSet(chrM), "chrM.fa")
}


# build index using base
if (!dir.exists("mm10_base_index")) {
  Sys.setenv(b2b="/Users/sr286/Dropbox/OSX_software/bowtie2-2.3.5.1-macos-x86_64/bowtie2-build")
  dir.create("mm10_base_index")
  Sys.setenv(base_fa=base_fa)
  system("$b2b --threads 5 -f $base_fa mm10_base_index/mm10_base")
}

# create kmer fasta file, I couldn't get bowtie2 -F flag to work!
if (!file.exists("chrM_kmers.fa")) {
  m <- readDNAStringSet("chrM.fa") %>% as.character()
  m_split <- sapply(seq(1, nchar(m)-34,1), function(i) substr(m, i, i+33))
  names(m_split) <- paste("kmer", 1:length(m_split), sep="")
  dna <- DNAStringSet(m_split)
  writeXStringSet(dna, "chrM_kmers.fa")
}

# align kmers of chrM
dir.create("chrM_homolog")
Sys.setenv(b2="/Users/sr286/Dropbox/OSX_software/bowtie2-2.3.5.1-macos-x86_64/bowtie2")
system("$b2 -x mm10_base_index/mm10_base -f -U chrM_kmers.fa -a -S chrM_homolog/out.sam")

Sys.setenv(st="/Users/sr286/Dropbox/OSX_software/samtools-1.3/samtools")
system("$st view -u -F 4  chrM_homolog/out.sam -o chrM_homolog/out.bam")

Sys.setenv(bt="/opt/local/bin/bedtools")
system("$bt bamtobed -i chrM_homolog/out.bam > chrM_blacklist.bed")

# call peaks
#Sys.setenv(mc2="/Users/sr286/opt/anaconda3/bin/macs2")
#system("$mc2 callpeak -t chrM_homolog/out.bam -f BAM --outdir chrM_homolog/ --nomodel --extsize 50 -g 1.87e9 -p .1")

system("cut -f 1,2,3 chrM_blacklist.bed > temp1.bed")
system("cut -f 1,2,3 mm10-blacklist.v2.bed > temp2.bed")
system("cat temp1.bed temp2.bed > temp.bed")
system("$bt sort -i temp.bed > joint_blacklist.bed")
system("rm temp.bed temp1.bed temp2.bed")

