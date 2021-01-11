workflow idr2isre {

    File cell_type_file = "cell_type.txt"
    File isre_bed = "../../all_ISRE_sorted.bed"

    Array[String] a = read_lines(cell_type_file)
    String cell_type = a[0]

    scatter (cell_type in a) {

      call moveidrFroms3 {
                  input:
                    cell_type=cell_type
      }

      call IntersectISRE {
                  input:
                   peaks_file=moveidrFroms3.bed_file,
                   all_ISRE_file=isre_bed
      }

      call moveISRETos3 {
              input:
                cell_type = cell_type,
                peaks_file=IntersectISRE.intersect_file
      }
    }
}

#############################

task moveidrFroms3 {
    String cell_type
    String s3_container = "yoshida-atacseq/idr/"

    command {
      aws s3 cp s3://yoshida-atacseq/idr/${cell_type}.bed out.bed
    }

    output {
        File bed_file = "out.bed"
    }
}

task IntersectISRE {
  File peaks_file
  File all_ISRE_file

  Int w = 10
  String ll="{"
  String rr="}"

  command {
    # here I assume that the ISRE file is already sorted!
    # windows based on summit +/- w
    #bedtools sort -i ${peaks_file} | cut -f 1,2,3,10 > C.bed
    # awk '${ll}print $1 "\t" $2+$4-${w} "\t" $2+$4+${w}${rr}' C.bed  > summit_windows.bed
    # OR use all the peak
    bedtools sort -i ${peaks_file} | cut -f 1,2,3 > summit_windows.bed

    bedtools sort -i summit_windows.bed > summit_windows_sorted.bed
    bedtools intersect -a ${all_ISRE_file} -b summit_windows_sorted.bed -wa -u  > intersect.bed
  }
  output {
    File intersect_file = "intersect.bed"
  }
}

task moveISRETos3 {
    String cell_type
    File peaks_file
    String s3_container = "yoshida-atacseq/ISRE/"

    command {
        aws s3 cp "${peaks_file}" "s3://${s3_container}${cell_type}_ISRE.bed"
    }
}
