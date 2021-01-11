workflow bam2peaks {

    # bam file is assumed to be in s3 with form accession.bam
    # accession file is a single column of accessions
    # groupName is used to create the S3 folder in which peak files are
    #   placed and the filenames within that folder.  Read in from file.

    # output files:  pooled peaks, replicate peaks
    # if there is a single accession, then replicate peaks are pseudo peaks

    File accession_file
    File group_name_file
    File genome

    Array[String] a = read_lines(accession_file)
    Array[String] g = read_lines(group_name_file)
    String groupName = g[0]

    # create narrowPeak for each accession
    scatter (accession in a) {

        call moveBamFroms3 {
                input:
                  accession=accession
        }

        call Bam2Tag {
                input:
                 accession=accession,
                 bam_file=moveBamFroms3.file
        }
    }

    # create pooled tag files, and if needed pseudo-replicates
    call ProcessTags {
      input:
        a=a,
        tag_files=Bam2Tag.tag_files,
    }

    Array[String] processed_accessions = flatten([["pooled"], ProcessTags.accessions])
    Array[File] tag_files = flatten([[ProcessTags.pooled], ProcessTags.replicates])

    scatter (i in range(length(tag_files))) {

        String current_accession = processed_accessions[i]

        call Tag2Peaks {
          input:
            accession=current_accession,
            tag_file=tag_files[i]
        }

        call FilterPeaks {
          input:
             accession=current_accession,
             macs2_files=Tag2Peaks.files,
             genome=genome
        }

        call movePeaksTos3 {
            input:
              groupName=groupName,
              accession=current_accession,
              filenames=FilterPeaks.files
        }
    }
}

#############################
task moveBamFroms3 {
    String accession
    String s3_container

    command {
        aws s3 cp "s3://${s3_container}${accession}.bam" "${accession}.bam"
    }
    output {
        File file="${accession}.bam"
    }
}

# convert bam file to tag.align
task Bam2Tag {

  String accession
  File bam_file

  String lcb="{"
  String rcb="}"

  command {
    # fixing read names to match across pairs
    samtools view -h ${bam_file} | sed '/^@/!s/^SRR[0-9]*\./SRR/' | sed '/^@/!s/\.[1-2]//' | samtools view -h - -o temp.bam
    # sorting indexing
    samtools sort -n temp.bam -o temp2.bam
    # creating bampe file
    bedtools bamtobed -bedpe -mate1 -i temp2.bam  > temp3.bedpe

    # convert bedpe to tagalign and shifts tn5 cuts (splits reads into separate lines)
    cat temp3.bedpe | awk 'BEGIN${lcb}OFS="\t"${rcb}${lcb}printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10${rcb}' > temp3.bedpe.temp3.tagalign

    cat temp3.bedpe.temp3.tagalign | awk -F $'\t' 'BEGIN ${lcb}OFS = FS${rcb}${lcb} if ($6 == "+") ${lcb}$2 = $2 + 4${rcb} else if ($6 == "-") ${lcb}$3 = $3- 5${rcb} print $0${rcb}' > ${accession}.tagalign
  }
  # problems with docker memory I think!  samtools sort wouldn't run
  #runtime {
  #    docker: 'sleviyang/genomics:ngs_utilities'
  #}
  output {
      File tag_files="${accession}.tagalign"
  }
}

# collect tag files, produce pooled tag file and pseudo replicates if only one tag file
task ProcessTags {
  Array[File] tag_files
  Array[String] a
  File pshuf_and_split

  String lcb="{"
  String rcb="}"

  Int ll = length(tag_files)

  command {
    if (( ${ll}==1 ))
    then
      # pooled file is the tag file and create pseudo replicates
      mv "${tag_files[0]}" "pooled.tagalign"

      cat "pooled.tagalign" | sed 's/\n/\t/' > joined.bedpe
      python3 ${pshuf_and_split} joined.bedpe pseudo1.bedpe pseudo2.bedpe

      # convert to tagalingn
      cat pseudo1.bedpe | awk 'BEGIN${lcb}OFS="\t"${rcb}${lcb}printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10${rcb}' \
        > pseudo1.tagalign
      cat pseudo2.bedpe | awk 'BEGIN${lcb}OFS="\t"${rcb}${lcb}printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10${rcb}' \
        > pseudo2.tagalign
    else
      cat ${sep=" " tag_files} > "pooled.tagalign"
      all_tag_files=(${sep=" " tag_files})
      touch tag_filenames.txt
      for i in $(seq 0 1 $(( ${ll}-1 )) )
      do
        ln $${lcb}all_tag_files[i]${rcb} "ProcessTagsTemp_$i"
        echo "ProcessTagsTemp_$i" >> tag_filenames.txt
      done
    fi
  }
  output {
    File pooled = "pooled.tagalign"
    Array[File] replicates = if (ll==1) then ["pseudo1.tagalign", "pseudo2.tagalign"] else read_lines("tag_filenames.txt")
    Array[String] accessions = if (ll==1) then ["${a[0]}_pseudo1", "${a[0]}_pseudo2"] else a
  }
}

# run macs2
task Tag2Peaks {
  String accession
  File tag_file

  String gs="1.87e9"
  # ENCODE VALUES
  String pval_thresh="0.10"
  String shiftsize="75"
  String extsize="150"

  String lcb="{"
  String rcb="}"

  command {
    # creates ${accession}_peaks.narrowPeak, _control_lambda.bdg, _treat_pileup.bdg
    macs2 callpeak -t ${tag_file} -f BED -n ${accession} -g ${gs} -p ${pval_thresh} \
    --shift ${shiftsize} --extsize ${extsize} --nomodel -B --SPMR --keep-dup all --call-summits

  }
  output {
    Array[File] files=["${accession}_peaks.narrowPeak"]
  }
}

# sort and remove blacklist
task FilterPeaks {
   String accession
   Array[File] macs2_files
   File genome
   File blacklist

   File narrowPeak=macs2_files[0]

   String lcb="{"
   String rcb="}"

   command {
     sort -k 8gr,8gr ${narrowPeak} | awk 'BEGIN${lcb}OFS="\t"${rcb}${lcb}$4="Peak_"NR ; print $0${rcb}' > narrowPeak_sorted.bed
     # remove blacklist from peak calls
     bedtools intersect -v -a narrowPeak_sorted.bed -b ${blacklist} > ${accession}_narrowPeak.bed

   }
   output {
     Array[File] files=["${accession}_narrowPeak.bed"]
   }

}

task movePeaksTos3 {
    String groupName
    String accession
    Array[File] filenames
    String s3_container

    command {
        aws s3 cp "${filenames[0]}" "s3://${s3_container}${groupName}/${accession}_narrowPeak.bed"
    }
}
