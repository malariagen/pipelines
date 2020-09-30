version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
## This initial version of the pipeline is designed to ONLY work on one sample
## It can take a list of input_crams, input_bams or input_fastqs (paired).
## If more than one of these lists of files are provided, the pipeline will use in order:
## input_crams first, input_bams second (if input_crams not provided), and input_fastqs lastly.
##

import "../../../structs/gcp/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/gcp/ReadBackedPhasingTasks.wdl" as Tasks

workflow ReadBackedPhasing {
  String pipeline_version = "0.0.0"

  input {
    File sample_id
    String output_basename
    File input_bam
    File input_bam_index
    File input_vcf
    Array[File] phased_sites_vcfs
    Array[File] phased_sites_indicies

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }
  # Step 1: Genotype data preparation
  call Tasks.SelectVariants {
    input:
      input_vcf = input_vcf,
      phased_sites_vcfs = phased_sites_vcfs,
      phased_sites_indicies = phased_sites_indicies,
      output_basename = output_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Step 2: WhatsHap phase
  call Tasks.WhatsHapPhase {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      subset_vcf = SelectVariants.subset_vcf,
      output_basename = output_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Step 3: WhatsHap stats
  call Tasks.WhatsHapStats {
    input:
      phased_vcf = WhatsHapPhase.phased_vcf,
      output_basename = output_basename,
      runTimeSettings = runTimeSettings
  }

  output {
   # TODO: determine outputs needed (stats etc.)
    File output_vcf = WhatsHapPhase.output_vcf
    File whats_hap_stats = "~{output_basename}.whatshap_stats"
    File subset_vcf = SelectVariants.subset_vcf
  }
}


