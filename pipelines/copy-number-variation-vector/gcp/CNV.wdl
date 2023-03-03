version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.


 import "../../structs/gcp/RunTimeSettings.wdl"
 import "../../structs/ReferenceSequence.wdl"
# import "../../../pipelines/phasing-vector/gcp/ReadBackedPhasing.wdl" as ReadBackedPhasing
# import "../../../pipelines/phasing-vector/gcp/StatisticalPhasing.wdl" as StatisticalPhasing

workflow CNV {
  String pipeline_version = "1.0.0"

  input {
    String project_id
    Array[String] sample_ids
    Array[File] input_bams
    Array[File] input_bam_indices
    Array[File] sample_zarrs
    Array[File] sample_vcfs = []
    Array[File] sample_vcf_indices = []
    File called_sites_zarr
    File phased_sites_zarr
    Array[String] chromosome_list
    Array[File] genetic_maps
    Array[File] interval_lists

    Array[File]? haplotype_reference_panels
    Array[File]? haplotype_reference_panel_indices

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
    String runtime_zones
  }

  # TODO: extract sample_id, sample_bam, and sample_vcf information from the maifest file (or inputs)

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

    # Array[File] output_vcf = StatisticalPhasing.output_vcf
    # Array[File] output_vcf_index = StatisticalPhasing.output_vcf_index
    # Array[File] output_zarrs = StatisticalPhasing.zarr_output

    # Array[Array[File]] read_back_phased_sample_vcfs = ReadBackedPhasing.phased_sample_vcf
    # Array[Array[File]] read_back_phased_sample_vcf_indices = ReadBackedPhasing.phased_sample_vcf_index
    # Array[Array[File]] whats_hap_stats_tsvs = ReadBackedPhasing.whats_hap_stats_tsv
    # Array[Array[File]] whats_hap_blocks_gtfs = ReadBackedPhasing.whats_hap_blocks_gtf
  }
}

task hello {
  command {
    echo "Hello, World!"
  }
  output {
    String message = read_string(stdout())
  }
}
