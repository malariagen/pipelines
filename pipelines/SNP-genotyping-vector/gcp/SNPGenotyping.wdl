version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Mosquito SNP Genotyping Pipeline as described in
## https://github.com/malariagen/pipelines/blob/991b0328ea8027bc6b1137f893a5340e27c8e87c/docs/specs/snp-genotyping-vector.md
## This is an initial proof of concept implementation.  It is designed to ONLY work on one sample, with ONLY its
## list of input fastqs.  It is currently implemented to run using Cromwell with a google cloud platform backend.
##

import "../../../structs/gcp/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/gcp/SNPGenotypingTasks.wdl" as Tasks

workflow SNPGenotyping {
  String pipeline_version = "1.0.1"

  input {
    String sample_id
    String output_basename = sample_id
    File input_bam
    File input_bam_index
    File alleles_vcf
    File alleles_vcf_index

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  call Tasks.UnifiedGenotyper {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      alleles_vcf = alleles_vcf,
      alleles_vcf_index = alleles_vcf_index,
      output_vcf_filename = output_basename + ".vcf",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.VcfToZarr {
    input:
      input_vcf = UnifiedGenotyper.output_vcf,
      input_vcf_index = UnifiedGenotyper.output_vcf_index,
      sample_id = sample_id,
      output_zarr_file_name = output_basename + ".zarr",
      output_log_file_name = output_basename + ".log",
      runTimeSettings = runTimeSettings
  }

  output {
    File output_vcf = UnifiedGenotyper.output_vcf
    File output_vcf_index = UnifiedGenotyper.output_vcf_index
    File zarr_output = VcfToZarr.zarr_output
  }
  meta {
    allowNestedInputs: true
  }
}
