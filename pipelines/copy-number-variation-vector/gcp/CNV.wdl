version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.


#  import "../../structs/gcp/RunTimeSettings.wdl"
#  import "../../structs/ReferenceSequence.wdl"
# import "../../../pipelines/phasing-vector/gcp/ReadBackedPhasing.wdl" as ReadBackedPhasing
# import "../../../pipelines/phasing-vector/gcp/StatisticalPhasing.wdl" as StatisticalPhasing

workflow CNV {
  String pipeline_version = "1.0.0"

  input {
    String project_id
    File sample_bam
  }

  

  # This is a wdl hack to create a pseudo None
  if (false) {
    File? none = "None"
  }

  meta {
    allowNestedInputs: true
  }
 
  call hello

  output {
    String test = hello.message
  }
}

task hello {
  command {
    echo "Hello, World!"
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1677557222"
    memory: "4 GiB"
    preemptible: 3
  }
  output {
    String message = read_string(stdout())
  }
}
