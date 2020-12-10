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
import "../../../tasks/gcp/Tasks.wdl" as Tasks
import "../../../tasks/gcp/ReadBackedPhasingTasks.wdl" as ReadBackedPhasingTasks

workflow ReadBackedPhasing {
  String pipeline_version = "0.0.0"

  input {
    File sample_id
    String output_basename
    File input_bam
    File input_bam_index
    File sample_zarr
    File called_sites_zarr
    File phased_sites_zarr
    String contig

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }
  # Step 1: Genotype data preparation
  call ReadBackedPhasingTasks.SelectVariants {
    input:
      sample_zarr = sample_zarr,
      called_sites_zarr = called_sites_zarr,
      phased_sites_zarr = phased_sites_zarr,
      output_basename = output_basename,
      contig = contig,
      runTimeSettings = runTimeSettings
  }
  # bgzip and index the phased_vcf (Needed by WhatsHap phase)
  call Tasks.BgzipAndTabix as BgzipAndTabixSelectVariantsVcf {
    input:
      input_vcf = SelectVariants.subset_vcf,
      output_basename = output_basename + ".subset",
      runTimeSettings = runTimeSettings
  }
  # Step 2: WhatsHap phase
  call ReadBackedPhasingTasks.WhatsHapPhase {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      subset_vcf = BgzipAndTabixSelectVariantsVcf.vcf,
      subset_vcf_index = BgzipAndTabixSelectVariantsVcf.vcf_index,
      output_filename = output_basename + ".phased.vcf.gz",
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # index the phased_vcf (needed for bcftools merge in statistical phasing pipeline)
  call Tasks.Tabix as TabixPhasedVcf {
    input:
      input_file = WhatsHapPhase.phased_vcf,
      runTimeSettings = runTimeSettings
  }
  # Step 3: WhatsHap stats
  call ReadBackedPhasingTasks.WhatsHapStats {
    input:
      phased_vcf = WhatsHapPhase.phased_vcf,
      output_basename = output_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    File subsetted_sample_vcf = SelectVariants.subset_vcf
    File phased_sample_vcf = WhatsHapPhase.phased_vcf
    File phased_sample_vcf_index = TabixPhasedVcf.index_file
    File whats_hap_stats_tsv = WhatsHapStats.whats_hap_stats_tsv
    File whats_hap_blocks_gtf = WhatsHapStats.whats_hap_blocks_gtf
  }
}


