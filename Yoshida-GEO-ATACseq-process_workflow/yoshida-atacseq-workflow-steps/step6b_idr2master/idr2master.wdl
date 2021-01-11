workflow idr2master {

    call formMaster

    call moveMasterTos3 {
            input:
              master_nonoverlapping_file=formMaster.master_nonoverlapping_file
    }
}

#############################

task formMaster {
    File filter_peaks_python = "filter_for_nonintersecting_idr.py"
    Int we = 250
    String ll="{"
    String rr="}"

    command {
      # download all the idr bed files
      mkdir temp_idr
      aws s3 cp s3://yoshida-atacseq/idr/ temp_idr/ --recursive

      # idr format (https://github.com/nboley/idr)
      # chr13	90922592	90923388	.	1000	.	24.30880	216.85700	211.80500	392	5.00	5.00
      # 90922592	90923504	19.12860	389	90922595	90923372	24.52260	392

      cat temp_idr/*  > master.bed
      awk '${ll}print $1 "\t"  $2+$10-${we} "\t" $2+$10+${we} "\t" $9 ${rr}' master.bed > master_window.bed
      sort -k1,1 -k2,2n master_window.bed > master_window_sorted.bed

      # run python script to filter for non-overlapping isre
      python3 ${filter_peaks_python} master_window_sorted.bed master_nonoverlap.bed
    }

    output {
        File master_nonoverlapping_file = "master_nonoverlap.bed"
    }
}


task moveMasterTos3 {
    File master_nonoverlapping_file
    String s3_container = "yoshida-atacseq/master/"

    command {
        aws s3 cp "${master_nonoverlapping_file}" "s3://${s3_container}master_idr_nonoverlapping.bed"

    }
}
