version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-LaneletAlignmentTask-vector.md
## This initial version of the pipeline is designed to ONLY work on one sample
## It can take a list of input_crams, input_bams or input_fastqs (paired).
## If more than one of these lists of files are provided, the pipeline will use in order:
## input_crams first, input_bams second (if input_crams not provided), and input_fastqs lastly.
##

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/ImportTasks.wdl" as ImportTask
import "../../../tasks/farm5/ShortReadAlignmentTasks.wdl" as SampleAlignmentTask
import "../../../pipelines/import-short-read-lanelet-alignment-vector/farm5/ImportShortReadLaneletAlignment.wdl" as LaneletAlignmentWorkflow

workflow ImportShortReadAlignment {
  String pipeline_version = "1.0.0"

  input {
    File batch_sample_manifest_file

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  call ImportTask.BatchSplitUpInputFile as BatchSplitUpInputFile {
    input:
      batch_sample_manifest_file = batch_sample_manifest_file,
      runTimeSettings = runTimeSettings
  }


  scatter (per_sample_manifest_file in BatchSplitUpInputFile.per_sample_manifest_files) {

    Array[Array[String]] lanelet_infos = read_tsv(per_sample_manifest_file)
    String sample_id = lanelet_infos[0][0]

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
        output_bam_basename = sample_id + ".fixmate.bam",
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
  }

  output {
    Array[File] indel_realigner_interval_list_file = RealignerTargetCreator.output_interval_list_file
    Array[File] output_bam = FixMateInformation.output_bam
    Array[File] output_bam_index = FixMateInformation.output_bam_index
    Array[File] validate_samfile_report_file = ValidateSamFile.report_file
    Array[File] samtools_stats_report_file = SamtoolsStats.report_file
    Array[File] samtools_idxstats_report_file = SamtoolsIdxStats.report_file
    Array[File] samtools_flagstat_report_file = SamtoolsFlagStat.report_file
    Array[File] callable_loci_summary_file = GatkCallableLoci.summary_file
  }
}
