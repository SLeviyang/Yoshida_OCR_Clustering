workflow peaks2idr {

    # cell_type_file should contain a single cell type.  The
    # cell_type is used to access the s3 directory
    # s3://yoshida-atacseq/peaks/$cell_type/

    File cell_type_file

    Array[String] a = read_lines(cell_type_file)
    String cell_type = a[0]

    call movePeaksFroms3 {
                input:
                  cell_type=cell_type
    }

    call ProcessPeaks {
                input:
                 pooled_file=movePeaksFroms3.pooled_file,
                 replicate_files=movePeaksFroms3.replicate_files
    }

    call movePeaksTos3 {
            input:
              cell_type=cell_type,
              peaks_file=ProcessPeaks.peaks_file
    }
}

#############################

task movePeaksFroms3 {
    String cell_type
    String s3_container

    command {
      mkdir pooled
      aws s3 cp s3://yoshida-atacseq/peaks/${cell_type}/ pooled/  --recursive --include "pooled*" --exclude "SRR*"
      cp pooled/* pooled.bed
      rm -r pooled/
      mkdir replicates
      aws s3 cp s3://yoshida-atacseq/peaks/${cell_type}/ replicates/  --recursive --include "SRR*" --exclude "pooled*"
    }

    output {
        File pooled_file = "pooled.bed"
        Array[File] replicate_files = glob("replicates/*")
    }
}

task ProcessPeaks {
  File pooled_file
  Array[File] replicate_files
  File idr

  Int nr = length(replicate_files)
  String ll="{"
  String rr="}"

  command {
    #${idr} --samples ${replicate_files[0]} ${replicate_files[1]}  \
    #    --peak-list ${pooled_file} \
    #    --output-file out.bed -i 0.05

    rfn=( ${sep=" " replicate_files} )
    indices="$( seq 0 1 $(( ${nr} - 1 )) )"

    # run idr on each pair and find pair with largest number of idr peaks
    npeaks_max=-1
    touch replicate_pair_peak_counts.txt
    for i in $indices
    do
      for j in $indices
      do
        if [[ $i < $j ]]
        then
          ${idr} --samples $${ll}rfn[i]${rr} $${ll}rfn[j]${rr}  \
              --peak-list ${pooled_file} \
              --output-file out_"r$i"_"r$j".bed -i 0.05
          npeaks=$(wc -l out_"r$i"_"r$j".bed | awk '${ll}print $1${rr}')
          echo -e "$i\t$j\t$npeaks" >> replicate_pair_peak_counts.txt
          if [[ $npeaks > $npeaks_max ]]
          then
            max_i=$i
            max_j=$j
            npeaks_max=$npeaks
          fi
        fi
      done
    done

    echo "$max_i $max_j $npeaks_max" > max.txt
    cp out_"r$max_i"_"r$max_j".bed final_peaks.bed
  }
  output {
    File peaks_file = "final_peaks.bed"
  }
}

task movePeaksTos3 {
    String cell_type
    File peaks_file
    String s3_container

    command {
        aws s3 cp "${peaks_file}" "s3://${s3_container}${cell_type}.bed"
    }
}
