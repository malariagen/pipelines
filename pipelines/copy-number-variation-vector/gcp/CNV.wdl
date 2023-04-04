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
    # windowed coverage inputs
    Array[File] input_bams
    Array[String] sample_names
    #String scripts_folder="/cnv/scripts"
    String output_dir="coverage"
    Int interval = 300
    Int window_size = 300
    Int min_qual = 10
    # coverage summary stats inputs
    Float accessibility_threshold = 0.9
    Float mapq_threshold = 0.5
    File accessibility_mask_file
    File mapq_file
    File sample_manifest
    File gc_content_file
    String sample_group_id
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
        input_bam = input_bams[idx],
        sample_name = sample_names[idx],
        output_dir = output_dir,
        interval = interval,
        window_size = window_size,
        min_qual = min_qual,
        accessibility_threshold = accessibility_threshold,
        mapq_threshold = mapq_threshold,
        accessibility_mask_file = accessibility_mask_file,
        mapq_file = mapq_file,
        sample_manifest = sample_manifest,
        gc_content_file = gc_content_file,
        sample_group_id = sample_group_id
    }
  }
  output {
    Array[File] hmm_outputs = HMM.output_gz
  }
}