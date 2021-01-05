version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements Mospquito Phasing Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/master/docs/specs/phasing-vector.md
##.
##

import "../../../structs/gcp/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../pipelines/phasing-vector/gcp/ReadBackedPhasing.wdl" as ReadBackedPhasing
import "../../../pipelines/phasing-vector/gcp/StatisticalPhasing.wdl" as StatisticalPhasing

workflow Phasing {
  String pipeline_version = "0.0.0"

  input {
    String project_id
    Array[String] sample_ids
    Array[File] input_bams
    Array[File] input_bam_indices
    Array[File] sample_zarrs
    Array[File] sample_vcfs = []
    File called_sites_zarr
    File phased_sites_zarr
    Array[String] chromosome_list
    Array[File] genetic_maps

    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # TODO: extract sample_id, sample_bam, and sample_vcf information from the maifest file (or inputs)

  # Step 1: Read-backed phasing

  # Scatter over chromosomes
  scatter(chr_idx in range(length(chromosome_list))) {
    String chromosome = chromosome_list[chr_idx]
    File genetic_map = genetic_maps[chr_idx]

    # Scatter over samples
    scatter(idx in range(length(sample_ids))) {
      # Run read-backed phasing on each sample (for each chromosome)
      call ReadBackedPhasing.ReadBackedPhasing as ReadBackedPhasing {
        input:
          sample_id = sample_ids[idx],
          output_basename = sample_ids[idx] + "_" + chromosome,
          input_bam = input_bams[idx],
          input_bam_index = input_bam_indices[idx],
          sample_zarr = sample_zarrs[idx],
#          sample_vcf = sample_vcfs[idx],
          called_sites_zarr = called_sites_zarr,
          phased_sites_zarr = phased_sites_zarr,
          contig = chromosome,
          reference = reference,
          runTimeSettings = runTimeSettings
      }
    }

    # combine samples
    Array[File] phased_sample_vcfs = ReadBackedPhasing.phased_sample_vcf
    Array[File] phased_sample_vcf_indices = ReadBackedPhasing.phased_sample_vcf_index

    # Step 2: Statistical phasing
    # run statistical phasing for all samples (for each chromosome)
    call StatisticalPhasing.StatisticalPhasing {
      input:
        project_id = project_id + "_" + chromosome,
        phased_sample_vcfs = phased_sample_vcfs,
        phased_sample_vcf_indices = phased_sample_vcf_indices,
        contig = chromosome,
        genetic_map = genetic_map,
        haplotype_reference_panel = haplotype_reference_panel,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  # Combine all chromosomes???

  output {
    Array[File] output_vcf = StatisticalPhasing.output_vcf
    Array[Array[File]] output_zarr_files = StatisticalPhasing.output_zarr_files

    Array[Array[File]] read_back_phased_sample_vcfs = ReadBackedPhasing.phased_sample_vcf
    Array[Array[File]] read_back_phased_sample_vcf_indices = ReadBackedPhasing.phased_sample_vcf_index
    Array[Array[File]] whats_hap_stats_tsvs = ReadBackedPhasing.whats_hap_stats_tsv
    Array[Array[File]] whats_hap_blocks_gtfs = ReadBackedPhasing.whats_hap_blocks_gtf
  }
}
