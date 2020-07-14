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

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/Alignment.wdl" as Alignment
import "../../../tasks/farm5/ShortReadAlignmentTasks.wdl" as Tasks

workflow ShortReadAlignment {
  String pipeline_version = "1.0.0"

  input {
    String sample_id
    String output_basename = sample_id
    File input_file

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  call Tasks.SplitUpInputFile as SplitUpInputFile {
    input:
      input_file = input_file,
      sample_id = sample_id,
      runTimeSettings = runTimeSettings
  }

  scatter (lanelet_file in SplitUpInputFile.lanelet_files) {

    Array[Array[String]] tmp = read_tsv(lanelet_file)
    String run_ena = tmp[0][1]
    String irods_path = tmp[0][2]
    String bam_path = tmp[0][3]
    String cram_path = tmp[0][4]
    String read1_path = tmp[0][5]
    String read2_path = tmp[0][6]

    if (cram_path != "") {
      # cram_path is defined, let's use it
      # but first, see if we need to transfer it over (FTP)
      String cram_filename_without_ftp_prefix = sub(cram_path, "ftp://", "")
      Boolean cram_is_an_ftp_path = (cram_filename_without_ftp_prefix != cram_path)
      if (cram_is_an_ftp_path) {
        call Tasks.Ftp as FtpCram {
          input:
            input_string = cram_path,
            runTimeSettings = runTimeSettings
        }
      }
    }
    if (cram_path == "" &&  bam_path != "") {
      # cram_path is NOT defined, but bam_path is, let's use that.
      # but first, see if we need to transfer it over (FTP)
      String bam_filename_without_ftp_prefix = sub(bam_path, "ftp://", "")
      Boolean bam_is_an_ftp_path = (bam_filename_without_ftp_prefix != bam_path)
      if (bam_is_an_ftp_path) {
        call Tasks.Ftp as FtpBam {
          input:
            input_string = bam_path,
            runTimeSettings = runTimeSettings
        }
      }
    }
    if (cram_path == "" &&  bam_path == "") {
      # cram_path and bam_path are both NOT defined, we'll use read1_path and read2_path (the fastqs)
      String fastq1_filename_without_ftp_prefix = sub(read1_path, "ftp://", "")
      Boolean fastq1_is_an_ftp_path = (fastq1_filename_without_ftp_prefix != read1_path)
      if (fastq1_is_an_ftp_path) {
        call Tasks.Ftp as FtpFastq1 {
          input:
            input_string = read1_path,
            runTimeSettings = runTimeSettings
        }
      }

      String fastq2_filename_without_ftp_prefix = sub(read2_path, "ftp://", "")
      Boolean fastq2_is_an_ftp_path = (fastq2_filename_without_ftp_prefix != read2_path)
      if (fastq2_is_an_ftp_path) {
        call Tasks.Ftp as FtpFastq2 {
          input:
            input_string = read2_path,
            runTimeSettings = runTimeSettings
        }
      }
    }

    call Alignment.Alignment as Alignment {
      input:
        sample_id = sample_id,
        read_group_id = run_ena,
        input_cram = select_first([FtpCram.output_file, cram_path]),
        input_bam = select_first([FtpBam.output_file, bam_path]),
        input_fastq1 = select_first([FtpFastq1.output_file, read1_path]),
        input_fastq2 = select_first([FtpFastq2.output_file, read2_path]),
        output_file_basename = run_ena,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  call Tasks.MergeSamFiles {
    input:
      input_files = Alignment.output_bam,
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
      known_indels_vcf = known_indels_vcf,
      output_interval_list_filename = output_basename + ".intervals",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.IndelRealigner {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      known_indels_vcf = known_indels_vcf,
      interval_list_file = RealignerTargetCreator.output_interval_list_file,
      output_bam_filename = output_basename + ".bam",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.FixMateInformation {
    input:
      input_file = IndelRealigner.output_bam,
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
    File output_bam_index = FixMateInformation.output_bam_index
    File validate_samfile_report_file = ValidateSamFile.report_file
    File samtools_stats_report_file = SamtoolsStats.report_file
    File samtools_idxstats_report_file = SamtoolsIdxStats.report_file
    File samtools_flagstat_report_file = SamtoolsFlagStat.report_file
    File callable_loci_summary_file = GatkCallableLoci.summary_file
  }
}
