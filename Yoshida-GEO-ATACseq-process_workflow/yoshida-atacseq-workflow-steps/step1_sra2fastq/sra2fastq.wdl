workflow sra2fastq {
    # accession file has a single column of accessions
    File accession_file
    String s3_container

    Array[String] a = read_lines(accession_file)

    scatter (ca in a) {

      call fastq_from_sra {
            input:
              accession=ca,
      }

      call move2s3 {
            input:
              fname=fastq_from_sra.out,
              s3_container=s3_container
      }
    }
}

task fastq_from_sra {
    String accession

    String modifier="-I --split-files"

    command {
        /sratoolkit/bin/fastq-dump ${modifier} "${accession}"
    }
    output {
          Array[File] out=["${accession}_1.fastq", "${accession}_2.fastq"]
    }
    runtime {
        docker: 'sleviyang/genomics:sratoolkit'
    }
}

task move2s3 {
    Array[File] fname
    String s3_container

    command {
         aws s3 mv ${fname[0]} "s3://${s3_container}"
         aws s3 mv ${fname[1]} "s3://${s3_container}"
    }
    output {
        String out = "s3://${s3_container}"
    }
}
