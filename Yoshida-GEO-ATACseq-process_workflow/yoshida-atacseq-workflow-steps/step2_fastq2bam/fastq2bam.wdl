workflow fastq2bam {
    # assumes that fastq files are on S3
    # accession file contains columns accession

    File accession_file
    String s3_container_fastq
    String s3_container_bam
    String genome
    Int threads_per_alignment # passed to bowtie2
    File picard_path
    File debug_bam_file

    Array[String] a = read_lines(accession_file)

    scatter (accession in a) {
        Array[String] fastq_filenames=["${accession}_1.fastq", "${accession}_2.fastq"]

        call movepairedFroms3 {
                input:
                  filename=fastq_filenames,
                  s3_container=s3_container_fastq
        }

        call alignpairs {
                input:
                 accession=accession,
                 fastq_files=movepairedFroms3.files,
                 genome=genome,
                 threads=threads_per_alignment
        }

        call FilterBam {
            input:
              accession=accession,
              picard_path=picard_path,
              sam_file=alignpairs.outfile,
              #sam_file=debug_bam_file,
              min_quality=30
        }


        call movepairedTos3 {
                input:
                  accession=accession,
                  filename=FilterBam.outfile,
                  s3_container=s3_container_bam
        }
    }
}

#############################

task movepairedFroms3 {
    Array[String] filename
    String s3_container

    command {
        aws s3 cp "s3://${s3_container}${filename[0]}" "${filename[0]}"
        aws s3 cp "s3://${s3_container}${filename[1]}" "${filename[1]}"
    }
    output {
        Array[File] files=["${filename[0]}", "${filename[1]}"]
    }
}

task alignpairs {
  String accession
  Array[File] fastq_files
  String genome
  Int threads

  command {
    # matches ENCODE ATAC-seq pipeline v1 except I use -X1000 instead of -X2000
    /bowtie2/bowtie2 -X1000 --mm -p ${threads} -x /bowtie2/genomes/${genome} -1 ${fastq_files[0]} -2 ${fastq_files[1]} -S "${accession}.sam"
  }
  runtime {
      docker: 'sleviyang/genomics:bowtie2'
  }
  output {
      File outfile="${accession}.sam"
  }
}

# assumes that picard.jar is available through picard_path and samtools is in PATH
task FilterBam {
  String accession
  File sam_file
  Int min_quality
  File picard_path

  File pp=picard_path

  command {
    samtools view -u -F 1804 -f 2 -q 30 "${sam_file}" -o temp1.bam
    java -jar "${pp}" SortSam -I temp1.bam  -O temp1_sorted.bam --SORT_ORDER queryname
    java -jar "${pp}" FixMateInformation  -I temp1_sorted.bam  -O temp1_sorted_fixmate.bam --ADD_MATE_CIGAR true --ASSUME_SORTED false
    java -jar "${pp}" MarkDuplicates -I temp1_sorted_fixmate.bam  -O temp1_sorted_fixmate_duprm.bam --REMOVE_DUPLICATES true -M "picard_metric_file.txt"

    samtools sort temp1_sorted_fixmate_duprm.bam -o "${accession}_full.bam"
    samtools index "${accession}_full.bam"

    # debug
    samtools view -h  -o "${accession}.bam" "${accession}_full.bam" chr1 chr2 \
    chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 \
    chr16 chr17 chr18 chr19 chrX
  }
  output {
    File outfile="${accession}.bam"
  }
}

task movepairedTos3 {
    String accession
    File filename
    String s3_container

    command {
        aws s3 cp "${filename}" "s3://${s3_container}${accession}.bam"
    }
}
