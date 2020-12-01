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

workflow StatisticalPhasing {
  String pipeline_version = "0.0.0"

  input {
    String project_id
    Array[File] sample_phased_vcfs
    Array[File] sample_phased_vcf_indices
#    File genetic_map

    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # Step 1: Merge VCFs
  call Tasks.MergeVcfs as MergeVcfs {
    input:
      sample_phased_vcfs = sample_phased_vcfs,
      sample_phased_vcf_indices = sample_phased_vcf_indices,
      project_id = project_id,
      runTimeSettings = runTimeSettings
  }

  # Step 2: ShapeIt4
  call Tasks.ShapeIt4 as ShapeIt4 {
    input:
      merged_vcf = MergeVcfs.merged_vcf,
      project_id = project_id,
      genetic_map = genetic_map,
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Possible Step: Ligate regions (?)
  call Tasks.LigateRegions {
    input:
      reference = reference,
      runTimeSettings = runTimeSettings
  }
  # Step: VCF to Zarr
    call Tasks.VcfToZarr {
      input:
        phased_vcf = ShapeIt4.phased_vcf
    }

  output {
  # TODO: determine outputs needed (stats etc.)
    File output_vcf = ShapeIt4.phased_vcf
    File zarr_output = VcfToZarr.zarr_output
  }
}
