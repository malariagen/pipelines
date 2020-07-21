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
import "../../../tasks/gcp/Alignment.wdl" as Alignment
import "../../../tasks/gcp/ShortReadAlignmentTasks.wdl" as Tasks

workflow ShortReadAlignment {
  String pipeline_version = "0.2.0"

  input {
    String sample_id
    Array[String] read_group_ids
    Array[File]? input_crams
    Array[File]? input_bams
    Array[File]? input_fastq1s
    Array[File]? input_fastq2s
    String output_basename

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  # Note that the order of input_crams, input_bams, and input_fastqs are all optional.
  # If more than one of them is provided, the order of precedence is:
  # crams first (if found)
  # bams next (if found)
  # fastqs last
  Array[File] input_files = select_first([input_crams, input_bams, input_fastq1s])

  scatter (idx in range(length(input_files))) {
    if (defined(input_crams)) {
      File input_cram = select_first([input_crams])[idx]
      call Alignment.Alignment as CramAlignment {
        input:
          sample_id = sample_id,
          read_group_id = read_group_ids[idx],
          input_cram = input_cram,
          output_file_basename = basename(input_cram, ".cram"),
          reference = reference,
          runTimeSettings = runTimeSettings
      }
    }
    if (!defined(input_crams) && defined(input_bams)) {
      File input_bam = select_first([input_bams])[idx]
      call Alignment.Alignment as BamAlignment {
        input:
          sample_id = sample_id,
          read_group_id = read_group_ids[idx],
          input_bam = input_bam,
          output_file_basename = basename(input_bam, ".bam"),
          reference = reference,
          runTimeSettings = runTimeSettings
      }
    }
    if (!defined(input_crams) && !defined(input_bams)) {
      File input_fastq1 = select_first([input_fastq1s])[idx]
      call Alignment.Alignment as FastqAlignment {
        input:
          sample_id = sample_id,
          read_group_id = read_group_ids[idx],
          input_fastq1 = input_fastq1,
          input_fastq2 = select_first([input_fastq2s])[idx],
          output_file_basename = basename(basename(input_fastq1, ".gz"), ".fastq"),
          reference = reference,
          runTimeSettings = runTimeSettings
      }
    }
    File alignment_output_bam = select_first([CramAlignment.output_bam, BamAlignment.output_bam, FastqAlignment.output_bam])
  }

  call Tasks.MergeSamFiles {
    input:
      input_files = alignment_output_bam,
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
      input_file = FixMateInformation.output_bam,
      report_filename = output_basename + ".validation_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.SamtoolsStats {
    input:
      input_file = FixMateInformation.output_bam,
      report_filename = output_basename + ".samtools_stats_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.SamtoolsIdxStats {
    input:
      input_bam = FixMateInformation.output_bam,
      input_bam_index = FixMateInformation.output_bam_index,
      report_filename = output_basename + ".samtools_idxstats_report.txt",
      runTimeSettings = runTimeSettings
  }

  call Tasks.SamtoolsFlagStat {
    input:
      input_bam = FixMateInformation.output_bam,
      report_filename = output_basename + ".samtools_flagstat_report.txt",
      runTimeSettings = runTimeSettings
  }

  call Tasks.GatkCallableLoci {
    input:
      input_bam = FixMateInformation.output_bam,
      input_bam_index = FixMateInformation.output_bam_index,
      summary_filename = output_basename + ".gatk_callable_loci.summary.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    File indel_realigner_interval_list_file = RealignerTargetCreator.output_interval_list_file
    File output_bam = FixMateInformation.output_bam
    File validate_samfile_report_file = ValidateSamFile.report_file
    File samtools_stats_report_file = SamtoolsStats.report_file
    File samtools_idxstats_report_file = SamtoolsIdxStats.report_file
    File samtools_flagstat_report_file = SamtoolsFlagStat.report_file
    File callable_loci_summary_file = GatkCallableLoci.summary_file
  }
}
