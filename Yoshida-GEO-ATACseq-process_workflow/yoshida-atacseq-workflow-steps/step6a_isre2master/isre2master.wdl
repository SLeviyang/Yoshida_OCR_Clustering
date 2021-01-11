workflow isre2master {

    # isre_bed contains all possible isres on mm10 genome (about 7 million)
    File isre_bed = "../../all_ISRE_sorted.bed"

    call formMaster

    call moveMasterTos3 {
            input:
              master_file=formMaster.master_file,
              master_nonoverlapping_file=formMaster.master_nonoverlapping_file
    }
}

#############################

task formMaster {
    File filter_peaks_python = "filter_for_nonintersecting_isre.py"

    command {
      # download all the isre bed files
      mkdir temp_ISRE
      aws s3 cp s3://yoshida-atacseq/ISRE/ temp_ISRE/ --recursive

      cat temp_ISRE/* | sort -k1,1 -k2,2n -k3,3n | uniq > master.bed

      # run python script to filter for non-overlapping isre
      python3 ${filter_peaks_python} master.bed master_nonoverlap.bed
    }

    output {
        File master_file = "master.bed"
        File master_nonoverlapping_file = "master_nonoverlap.bed"
    }
}


task moveMasterTos3 {
    File master_file
    File master_nonoverlapping_file
    String s3_container = "yoshida-atacseq/master/"

    command {
        aws s3 cp "${master_file}" "s3://${s3_container}master_isre.bed"
        aws s3 cp "${master_nonoverlapping_file}" "s3://${s3_container}master_isre_nonoverlapping.bed"

    }
}
