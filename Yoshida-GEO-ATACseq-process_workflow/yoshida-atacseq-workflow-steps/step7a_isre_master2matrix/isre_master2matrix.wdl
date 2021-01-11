workflow isre_master2cell_type_master {

    # cell_type_file contains the cell type
    #File cell_type_file = "cell_type.txt"

    #Array[String] g = read_lines(cell_type_file)
    #String cell_type = g[0]

    call loadMasterISRE
    call toMatrix {
      input:
        master_bed = loadMasterISRE.master_bed
    }

    call tos3 {
      input:
        matrix_file = toMatrix.matrix_file
    }



}

#############################
task loadMasterISRE {

    command {
      aws s3 cp s3://yoshida-atacseq/master/master_isre_nonoverlapping.bed master.bed
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
      mkdir temp_ISRE
      aws s3 cp s3://yoshida-atacseq/ISRE/ temp_ISRE/ --recursive

      mkdir zo
      all_files=( $(ls temp_ISRE) )
      for cf in $${ll}all_files[@]${rr}
      do
        cp ${master_bed} master.bed
        cp temp_ISRE/$cf current.bed
        bedtools closest -a ${master_bed} -b temp_ISRE/$cf -d -t first | \
          awk '${ll}if ($17==0) print "1"; else print "0"${rr}' > zo/$cf
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
    aws s3 cp ${matrix_file} s3://yoshida-atacseq/master/master_matrix.txt
  }
}
