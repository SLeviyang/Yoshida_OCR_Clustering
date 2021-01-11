workflow iidr_master2cell_type_master {

    call loadMaster
    call toMatrix {
      input:
        master_bed = loadMaster.master_bed
    }

    call tos3 {
      input:
        matrix_file = toMatrix.matrix_file
    }



}

#############################
task loadMaster {

    command {
      aws s3 cp s3://yoshida-atacseq/master/master_idr_nonoverlapping.bed master.bed
    }
    output {
      File master_bed = "master.bed"
    }
}

task toMatrix {
    File master_bed
    String ll = "{"
    String rr = "}"

    command {
      mkdir temp_idr
      aws s3 cp s3://yoshida-atacseq/idr/ temp_idr/ --recursive

      mkdir zo
      all_files=( $(ls temp_idr) )
      for cf in $${ll}all_files[@]${rr}
      do
        cp temp_idr/$cf current.bed
        cp ${master_bed} master.bed
        # NEED TO CHECK WINDOW WIDTHS, SOMETHING OFF HERE
        # IDR CONTAINS SAME PEAK SEVERAL TIMES!
        sort -k1,1 -k2,2n -k3,3n current.bed > current_sorted.bed
        bedtools closest -a ${master_bed} -b current_sorted.bed -d -t first | \
          awk '${ll}if ($25==0) print "1"; else print "0"${rr}' > zo/$cf
      done

      echo $${ll}all_files[@]${rr} > matrix.txt
      paste -d " " zo/* >> matrix.txt
    }
    output {
      File matrix_file = "matrix.txt"
    }
}

task tos3 {
  File matrix_file
  command {
    aws s3 cp ${matrix_file} s3://yoshida-atacseq/master/master_idr_matrix.txt
  }
}
