version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Statistical component of the Mospquito Phasing Pipeline as described in
## https://github.com/malariagen/pipelines/blob/master/docs/specs/phasing-vector.md
##

import "../../../structs/gcp/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/gcp/Tasks.wdl" as Tasks
import "../../../tasks/gcp/StatisticalPhasingTasks.wdl" as StatisticalPhasingTasks

workflow StatisticalPhasing {
  String pipeline_version = "1.0.0"

  input {
    String project_id
    Array[File] phased_sample_vcfs
    Array[File] phased_sample_vcf_indices
    String contig
    File genetic_map
    File interval_list

    #TODO - plug in haplotype_reference_panel
    File? haplotype_reference_panel

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # Step 1: Merge VCFs
  call StatisticalPhasingTasks.MergeVcfs as MergeVcfs {
    input:
      phased_sample_vcfs = phased_sample_vcfs,
      phased_sample_vcf_indices = phased_sample_vcf_indices,
      project_id = project_id,
      runTimeSettings = runTimeSettings
  }

  # Step 2: bgzip the merged VCF
  call Tasks.BgzipAndTabix {
    input:
      input_vcf = MergeVcfs.merged_vcf,
      output_basename = project_id + "_merged",
      runTimeSettings = runTimeSettings
  }

  # Step 3: ShapeIt4
  scatter(region in read_lines(interval_list)) {
    call StatisticalPhasingTasks.ShapeIt4 as ShapeIt4 {
      input:
        merged_vcf = BgzipAndTabix.vcf,
        merged_vcf_index = BgzipAndTabix.vcf_index,
        project_id = project_id,
        region = region,
        genetic_map = genetic_map,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
    call Tasks.Tabix as Tabix {
      input:
        input_file = ShapeIt4.region_phased_vcf,
        runTimeSettings = runTimeSettings
    }
  }

  # Step 4: Ligate regions
  call StatisticalPhasingTasks.LigateRegions as LigateRegions {
    input:
      region_phased_vcfs = Tabix.output_file,
      region_phased_vcfs_indicies = Tabix.output_file,
      interval_list = interval_list,
      project_id = project_id,
      runTimeSettings = runTimeSettings
  }

  # Step 5: Cohort VCF to Zarr
  call StatisticalPhasingTasks.CohortVcfToZarr {
    input:
      input_vcf = LigateRegions.phased_vcf,
      contig = contig,
      output_zarr_file_name = project_id + ".zarr",
      output_log_file_name = project_id + ".log",
      runTimeSettings = runTimeSettings
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File output_vcf = LigateRegions.phased_vcf
    File zarr_output = CohortVcfToZarr.zarr_output
  }
}
