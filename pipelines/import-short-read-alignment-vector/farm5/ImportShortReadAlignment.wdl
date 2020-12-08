version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-LaneletAlignmentTask-vector.md
## Imports multiple lanelet bams/crams for the same sample from IRODS.
## BWA aligns the lanelet bams/crams.  Processes each lanelet in parallel.
## Then merges the lanelets together into a single per-sample bam
## and IndelRealigns the bam and collects alignment stats.
## Lanelet refers a lane of a sample that has been multiplexed into multiple
## sequencing lanes.
## Expects that the per_sample_manifest_file contains only the lanelets for a single sample
## and each lanelet has its own row.
## Columns:  sample, run_ena, irods_path
##

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/ShortReadAlignmentTasks.wdl" as SampleAlignmentTask
import "../../../pipelines/import-short-read-lanelet-alignment-vector/farm5/ImportShortReadLaneletAlignment.wdl" as LaneletAlignmentWorkflow

workflow ImportShortReadAlignment {
  String pipeline_version = "1.0.0"

  input {
    String sample_id
    String output_basename
    File per_sample_manifest_file

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  call LaneletAlignmentWorkflow.ImportShortReadLaneletAlignment {
    input:
      sample_id = sample_id,
      per_sample_manifest_file = per_sample_manifest_file,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.MergeSamFiles {
    input:
      input_files = ImportShortReadLaneletAlignment.output_bam,
      output_filename = sample_id + ".mergesam.bam",
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.MarkDuplicates {
    input:
      input_bam = MergeSamFiles.output_file,
      output_filename = sample_id + ".markdup.bam",
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.RealignerTargetCreator {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      known_indels_vcf = known_indels_vcf,
      output_interval_list_filename = sample_id + ".intervals",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.IndelRealigner {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      known_indels_vcf = known_indels_vcf,
      interval_list_file = RealignerTargetCreator.output_interval_list_file,
      output_bam_filename = sample_id + ".indelrealign.bam",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.FixMateInformation {
    input:
      input_file = IndelRealigner.output_bam,
      output_bam_basename = sample_id + ".fixmate",
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.ValidateSamFile {
    input:
      input_file = FixMateInformation.output_bam,
      report_filename = sample_id + ".validation_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.SamtoolsStats {
    input:
      input_file = FixMateInformation.output_bam,
      report_filename = sample_id + ".samtools_stats_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.SamtoolsIdxStats {
    input:
      input_bam = FixMateInformation.output_bam,
      input_bam_index = FixMateInformation.output_bam_index,
      report_filename = sample_id + ".samtools_idxstats_report.txt",
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.SamtoolsFlagStat {
    input:
      input_bam = FixMateInformation.output_bam,
      report_filename = sample_id + ".samtools_flagstat_report.txt",
      runTimeSettings = runTimeSettings
  }

  call SampleAlignmentTask.GatkCallableLoci {
    input:
      input_bam = FixMateInformation.output_bam,
      input_bam_index = FixMateInformation.output_bam_index,
      summary_filename = sample_id + ".gatk_callable_loci.summary.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }


  output {
    File indel_realigner_interval_list_file = RealignerTargetCreator.output_interval_list_file
    File output_bam = FixMateInformation.output_bam
    File output_bam_index = FixMateInformation.output_bam_index
    File validate_samfile_report_file = ValidateSamFile.report_file
    File samtools_stats_report_file = SamtoolsStats.report_file
    File samtools_idxstats_report_file = SamtoolsIdxStats.report_file
    File samtools_flagstat_report_file = SamtoolsFlagStat.report_file
    File callable_loci_summary_file = GatkCallableLoci.summary_file
  }

  meta {
    allowNestedInputs: true
  }
}
