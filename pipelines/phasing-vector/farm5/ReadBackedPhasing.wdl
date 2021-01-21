version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2021
##
## This WDL pipeline implements the Read-Backed component of the Mospquito Phasing Pipeline as described in
## https://github.com/malariagen/pipelines/blob/master/docs/specs/phasing-vector.md
##

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/Tasks.wdl" as Tasks
import "../../../tasks/farm5/ShortReadAlignmentTasks.wdl" as ShortReadAlignmentTasks
import "../../../tasks/farm5/SNPGenotypingTasks.wdl" as SNPGenotypingTasks
import "../../../tasks/farm5/ReadBackedPhasingTasks.wdl" as ReadBackedPhasingTasks

workflow ReadBackedPhasing {
  String pipeline_version = "1.0.0"

  input {
    String sample_id
    String output_basename
    File input_bam
    File? input_bam_index
    File? sample_zarr
    File? sample_vcf
    File called_sites_zarr
    File phased_sites_zarr
    String contig

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  if (!defined(sample_zarr)) {
    # TODO - some sort of error if VCF not provided!
    call SNPGenotypingTasks.VcfToZarr {
      input:
        input_vcf = select_first([sample_vcf]),
        sample_id = sample_id,
        output_zarr_file_name = output_basename + ".zarr",
        output_log_file_name = output_basename + ".log",
        runTimeSettings = runTimeSettings
    }
  }

  # Step 1: Genotype data preparation
  call ReadBackedPhasingTasks.SelectVariants {
    input:
      sample_zarr = select_first([sample_zarr, VcfToZarr.zarr_output]),
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

  # If the bam index is not provided, generate it.
  if (!defined(input_bam_index)) {
    call ShortReadAlignmentTasks.SamtoolsIndex {
      input:
        input_file = input_bam,
        runTimeSettings = runTimeSettings
    }
  }

  # Step 2: WhatsHap phase
  call ReadBackedPhasingTasks.WhatsHapPhase {
    input:
      input_bam = select_first([SamtoolsIndex.output_file, input_bam]),
      input_bam_index = select_first([SamtoolsIndex.output_index_file, input_bam_index]),
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
      phased_vcf = TabixPhasedVcf.output_file,
      phased_vcf_index = TabixPhasedVcf.output_index_file,
      output_basename = output_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    File subsetted_sample_vcf = SelectVariants.subset_vcf
    File phased_sample_vcf = TabixPhasedVcf.output_file
    File phased_sample_vcf_index = TabixPhasedVcf.output_index_file
    File whats_hap_stats_tsv = WhatsHapStats.whats_hap_stats_tsv
    File whats_hap_blocks_gtf = WhatsHapStats.whats_hap_blocks_gtf
  }
}


