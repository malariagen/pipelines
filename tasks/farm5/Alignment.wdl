version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the lanelet-level Alignment workflow - a subworkflow of the Short Read Alignment Pipeline
## as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
##

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"
import "../../tasks/farm5/ShortReadAlignmentTasks.wdl" as Tasks

workflow Alignment {

  input {
    String sample_id
    String read_group_id
    File? input_cram
    File? input_bam
    File? input_fastq1
    File? input_fastq2
    String output_file_basename

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  if (defined(input_cram)) {
    call Tasks.CramToBam {
      input:
        input_file = select_first([input_cram]),
        output_filename = output_file_basename + ".bam",
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  if (defined(input_cram) || defined(input_bam)) {
    call Tasks.RevertSam {
      input:
        input_file = select_first([CramToBam.output_file, input_bam]),
        output_filename = output_file_basename + ".unmapped.bam",
        runTimeSettings = runTimeSettings
    }

    call Tasks.SamToFastq {
      input:
        input_file = RevertSam.output_file,
        output_fastq1_filename = output_file_basename + "_1.fastq",
        output_fastq2_filename = output_file_basename + "_2.fastq",
        runTimeSettings = runTimeSettings
    }
  }

  call Tasks.ReadAlignment {
    input:
      sample_id = sample_id,
      read_group_id = read_group_id,
      fastq1 = select_first([SamToFastq.output_fastq1, input_fastq1]),
      fastq2 = select_first([SamToFastq.output_fastq2, input_fastq2]),
      output_sam_basename = output_file_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call Tasks.ReadAlignmentPostProcessing {
    input:
      input_sam = ReadAlignment.output_sam,
      output_bam_basename = output_file_basename,
      runTimeSettings = runTimeSettings
  }

  call Tasks.SetNmMdAndUqTags {
    input:
      input_bam = ReadAlignmentPostProcessing.output_bam,
      output_bam_basename = output_file_basename,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    File output_bam = SetNmMdAndUqTags.output_bam
  }
}