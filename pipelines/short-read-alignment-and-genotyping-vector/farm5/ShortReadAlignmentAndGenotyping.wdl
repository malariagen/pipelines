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
import "../../../pipelines/short-read-alignment-vector/farm5/ShortReadAlignment.wdl" as ShortReadAlignment
import "../../../pipelines/SNP-genotyping-vector/farm5/SNPGenotyping.wdl" as SNPGenotyping

workflow ShortReadAlignmentAndGenotyping {
  String pipeline_version = "1.0.0"

  input {
    String sample_id
    String output_basename = sample_id
    File input_file

    File? known_indels_vcf
    File alleles_vcf
    File alleles_vcf_index

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  call ShortReadAlignment.ShortReadAlignment {
    input:
      sample_id = sample_id,
      output_basename = output_basename,
      input_file = input_file,
      known_indels_vcf = known_indels_vcf,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SNPGenotyping.SNPGenotyping {
    input:
      sample_id = sample_id,
      output_basename = output_basename,
      input_bam = ShortReadAlignment.output_bam,
      input_bam_index = ShortReadAlignment.output_bam_index,
      alleles_vcf = alleles_vcf,
      alleles_vcf_index = alleles_vcf_index,
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    File indel_realigner_interval_list_file = ShortReadAlignment.indel_realigner_interval_list_file
    File output_bam = ShortReadAlignment.output_bam
    File output_bam_index = ShortReadAlignment.output_bam_index
    File validate_samfile_report_file = ShortReadAlignment.validate_samfile_report_file
    File samtools_stats_report_file = ShortReadAlignment.samtools_stats_report_file
    File samtools_idxstats_report_file = ShortReadAlignment.samtools_idxstats_report_file
    File samtools_flagstat_report_file = ShortReadAlignment.samtools_flagstat_report_file
    File callable_loci_summary_file = ShortReadAlignment.callable_loci_summary_file
    File output_vcf = SNPGenotyping.output_vcf
    File output_vcf_index = SNPGenotyping.output_vcf_index
    File zarr_output = SNPGenotyping.zarr_output
  }
}
