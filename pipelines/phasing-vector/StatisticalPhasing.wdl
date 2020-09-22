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
import "../../../tasks/gcp/StatisticalPhasingTasks.wdl" as Tasks

workflow ReadBackedPhasing {
  String pipeline_version = "0.0.0"

  input {
    File sample_manifest
    Array[File] sample_bams
    Array[File] sample_vcfs
    File alleles_vcf
    File genetic_map

    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }
  # Step: Merge VCFs
  call Tasks.MergeVcfs {
    input:
      sample_manifest = sample_manifest,
      sample_bams = sample_bams,
      sample_vcfs = sample_vcfs,
      alleles_vcf = alleles_vcf,
      genetic_map = genetic_map,
      haplotype_reference_panel = haplotype_reference_panel,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Step: ShapeIt4
  call Tasks.ShapeIt4 {
    input:
      sample_manifest = sample_manifest,
      sample_bams = sample_bams,
      sample_vcfs = sample_vcfs,
      alleles_vcf = alleles_vcf,
      genetic_map = genetic_map,
      haplotype_reference_panel = haplotype_reference_panel,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Possible Step: Ligate regions (?)
  call Tasks.LigateRegions {
    input:
      sample_manifest = sample_manifest,
      sample_bams = sample_bams,
      sample_vcfs = sample_vcfs,
      alleles_vcf = alleles_vcf,
      genetic_map = genetic_map,
      haplotype_reference_panel = haplotype_reference_panel,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Step: VCF to Zarr
    call Tasks.VcfToZarr {
      input:
        sample_manifest = sample_manifest,
        sample_bams = sample_bams,
        sample_vcfs = sample_vcfs,
        alleles_vcf = alleles_vcf,
        genetic_map = genetic_map,
        haplotype_reference_panel = haplotype_reference_panel,
        reference = reference,
        runTimeSettings = runTimeSettings
    }

  output {
  # TODO: determine outputs needed (stats etc.)
    File output_vcf = ShapeIt4.output_vcf
    File zarr_output = VcfToZarr.zarr_output
  }
}
