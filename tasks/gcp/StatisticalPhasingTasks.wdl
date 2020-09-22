version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {

  }
  command {

  }
  runtime {

  }
  output {

  }
}


task ShapeIt4 {
# guts are not updated
  input {
    File input_file
    String sample_id
    RunTimeSettings runTimeSettings
  }

  String lanelet_file_prefix = "lanelet_"

  command {
    mkdir lanelet_temp
    cd lanelet_temp

    # Verify that the input_file begins with a header
    # should look like: `sample_id	run_ena	irods_path	bam_path	cram_path	read1_path	read2_path`
    head -1 ~{input_file} | grep '^sample_id\trun_ena'
    exitCode=$?
    if [ $exitCode != 0 ]; then
      echo "Input file ~{input_file} appears malformed"
      exit $exitCode
    fi

    # splits list of mappings into single files.  One line each.
    grep '^~{sample_id}\t' ~{input_file} | split -l 1 - ~{lanelet_file_prefix}
  }

  runtime {
    docker: runTimeSettings.lftp_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
      Array[File] lanelet_files = glob("lanelet_temp/~{lanelet_file_prefix}*")
  }
}

task LigateRegions {
  input {

  }
  command {

  }
  runtime {

  }
  output {

  }
}

task VcfToZarr {
# TODO: can we reuse the vcf to zarr from the genotyping pipeline?
  input {

  }
  command {

  }
  runtime {

  }
  output {

  }
}