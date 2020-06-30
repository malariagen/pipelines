version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
## This is an initial proof of concept implementation.  It is designed to ONLY work on one sample, with ONLY its
## list of input fastqs.  It is currently implemented to run using Cromwell with a google cloud platform backend.
##

import "../../structs/<PLATFORM>/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"
import "../../tasks/<PLATFORM>/ShortReadAlignmentTasks.wdl" as Tasks

workflow ShortReadAlignment {
  String pipeline_version = "0.1.0"

  input {
    String sample_id
    Array[String] read_group_ids
    Array[File] input_fastq1s
    Array[File] input_fastq2s
    String output_basename

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  scatter (idx in range(length(input_fastq1s))) {
    call Tasks.ReadAlignment {
      input:
        sample_id = sample_id,
        read_group_id = read_group_ids[idx],
        fastq1 = input_fastq1s[idx],
        fastq2 = input_fastq2s[idx],
        output_sam_basename = output_basename,
        reference = reference,
        runTimeSettings = runTimeSettings
    }

    call Tasks.ReadAlignmentPostProcessing {
      input:
        input_sam = ReadAlignment.output_sam,
        output_bam_basename = output_basename,
        runTimeSettings = runTimeSettings
    }

    call Tasks.SetNmMdAndUqTags {
      input:
        input_bam = ReadAlignmentPostProcessing.output_bam,
        output_bam_basename = output_basename,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  call Tasks.MergeSamFiles {
    input:
      input_files = SetNmMdAndUqTags.output_bam,
      output_filename = output_basename + ".bam",
      runTimeSettings = runTimeSettings
  }

  call Tasks.MarkDuplicates {
    input:
      input_bam = MergeSamFiles.output_file,
      output_filename = output_basename + ".bam",
      runTimeSettings = runTimeSettings
  }

  call Tasks.RealignerTargetCreator {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      output_interval_list_filename = output_basename + ".intervals",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.IndelRealigner {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      interval_list_file = RealignerTargetCreator.output_interval_list_file,
      output_filename = output_basename + ".bam",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.FixMateInformation {
    input:
      input_file = IndelRealigner.output_file,
      output_bam_basename = output_basename,
      runTimeSettings = runTimeSettings
  }

  call Tasks.ValidateSamFile {
    input:
      input_file = FixMateInformation.output_file,
      report_filename = output_basename + ".validation_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.SamtoolsStats {
    input:
      input_file = FixMateInformation.output_file,
      report_filename = output_basename + ".samtools_stats_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.GatkCallableLoci {
    input:
      input_bam = FixMateInformation.output_file,
      input_bam_index = FixMateInformation.output_index_file,
      summary_filename = output_basename + ".gatk_callable_loci.summary.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    Array[File] output_sams = ReadAlignment.output_sam
    Array[File] output_bams = ReadAlignmentPostProcessing.output_bam
    File indel_realigner_interval_list_file = RealignerTargetCreator.output_interval_list_file
    File output_bam = IndelRealigner.output_file
    File validate_samfile_report_file = ValidateSamFile.report_file
    File samtools_stats_report_file = SamtoolsStats.report_file
    File callable_loci_summary_file = GatkCallableLoci.summary_file
  }
}
