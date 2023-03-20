version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.


import "HMM.wdl" as HMM
workflow CNV {
  meta {
    description: "This is a pipeline for calling Copy Number Variants (CNVs) for a cohort of multiple samples. This inclludes an HMM step followed by coverage calls and target regions pipelines to improve acuracy."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    String project_id
    Array[File] input_bams
    String scripts_folder="/cnv/scripts"
  }

  

  # This is a wdl hack to create a pseudo None
  if (false) {
    File? none = "None"
  }

  # Scatter over samples
  scatter(idx in range(length(input_bams))) {
    # Run windowed coverage on each sample (for each chromosome)
    call HMM.HMM as HMM {
      input:
        input_bam = input_bams[idx]
    }
  }

  

  output {
    String test = ""
  }
}

# task hello {
#   command {
#     echo "Hello, World!"
#   }
#   runtime {
#     docker: "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1677557222"
#     memory: "4 GiB"
#     preemptible: 3
#   }
#   output {
#     String message = read_string(stdout())
#   }
# }

