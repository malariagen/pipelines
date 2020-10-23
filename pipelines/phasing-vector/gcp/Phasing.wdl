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
    Array[File] input_bams
    Array[File] input_bam_indicies
    Array[File] sample_zarrs
    File called_sites_zarr
    File phased_sites_zarr
    Array[String] chromosome_list
    File genetic_map # recombination rates

    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # TODO: extract sample_id, sample_bam, and sample_vcf information from the maifest file (or inputs)
  Array[String] sample_ids = []


  # Step 1: Read-backed phasing

  # Scatter over chormosomes
  scatter(chromosome in chromosome_list) {
    # Scatter over samples
    scatter(idx in range(length(sample_ids))) {
      # Run read-backed phasing on each sample (for each chromosome)
      call ReadBackedPhasing.ReadBackedPhasing as ReadBackedPhasing {
        input:
          sample_id = sample_ids[idx],
          input_bam = input_bams[idx],
          input_bam_index = input_bam_indicies[idx],
          sample_zarr = sample_zarrs[idx],
          called_sites_zarr = called_sites_zarr,
          phased_sites_zarr = phased_sites_zarr,
          contig = chromosome,
          reference = reference,
          runTimeSettings = runTimeSettings
      }
    }

    # combine samples
    Array[File] sample_phased_vcfs = ReadBackedPhasing.phased_vcf

    # Step 2: Statistical phasing
    # run statistical phasing for all samples (for each chromosome)
    call StatisticalPhasing.StatisticalPhasing {
      input:
      project_id = project_id,
      sample_phased_vcfs = sample_phased_vcfs,
      genetic_map = genetic_map,
      haplotype_reference_panel = haplotype_reference_panel,
      reference = reference,
      runTimeSettings = runTimeSettings
    }
  }

  # Combine all chromosomes

  output {
  # TODO: determine addtional outputs
    File output_vcf = StatisticalPhasing.output_vcf
    File zarr_output = StatisticalPhasing.zarr_output
  }
}
