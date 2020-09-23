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
import "../../../pipelines/phasing-vector/ReadBackedPhasing.wdl" as ReadBackedPhasing
import "../../../pipelines/pipelines/phasing-vector/StatisticalPhasing.wdl" as StatisticalPhasing

workflow Phasing {
  String pipeline_version = "0.0.0"

  input {
    String project_id
    File sample_manifest
    Array[File] sample_bams
    Array[File] sample_vcfs
    File alleles_vcf
    File genetic_map # recombination rates

    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # TODO: extract sample_id, sample_bam, and sample_vcf information from the maifest file (or inputs)
  Array[String] sample_ids = []

  # Step 1: Read-backed phasing
  # TODO: scatter over samples
  scatter(idx in range(length(sample_bams))) {
    call ReadBackedPhasing.ReadBackedPhasing as ReadBackedPhasing {
      input:
        sample_id = sample_ids[idx],
        sample_bam = sample_bams[idx],
        sample_vcf = sample_vcfs[idx],
        alleles_vcf = alleles_vcf,
        genetic_map = genetic_map, # probably not needed for this sub-pipeline
        haplotype_reference_panel = haplotype_reference_panel, # probably not needed for this sub-pipeline
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  Array[File] sample_phased_vcfs = ReadBackedPhasing.phased_vcf

  # Step 2: Statistical phasing
  call StatisticalPhasing.StatisticalPhasing {
    input:
    project_id = project_id,
    sample_phased_vcfs = sample_phased_vcfs,
    genetic_map = genetic_map,
    haplotype_reference_panel = haplotype_reference_panel,
    reference = reference,
    runTimeSettings = runTimeSettings
  }

  output {
  # TODO: determine addtional outputs
    File output_vcf = StatisticalPhasing.output_vcf
    File zarr_output = StatisticalPhasing.zarr_output
  }
}
