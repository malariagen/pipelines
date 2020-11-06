version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
## This version of the pipeline is designed to work on a batch of samples
## specified in a manifest with fields:
## sample_id, run_ena, irods_path
## It can take a list of input_crams, input_bams or input_fastqs (paired).
## If more than one of these lists of files are provided, the pipeline will use in order:
## input_crams first, input_bams second (if input_crams not provided), and input_fastqs lastly.
##

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/ImportTasks.wdl" as ImportTask
import "../../../pipelines/import-short-read-alignment-vector/farm5/ImportShortReadAlignment.wdl" as ImportShortReadAlignment
import "../../../pipelines/SNP-genotyping-vector/farm5/SNPGenotyping.wdl" as SNPGenotyping

workflow BatchImportShortReadAlignmentAndGenotyping {
  String pipeline_version = "1.0.0"

  input {
    File batch_sample_manifest_file

    File? known_indels_vcf
    File alleles_vcf
    File alleles_vcf_index

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
    String output_basename = sample_id

    call ImportShortReadAlignment.ImportShortReadAlignment {
      input:
        sample_id = sample_id,
        output_basename = output_basename,
        per_sample_manifest_file = per_sample_manifest_file,
        known_indels_vcf = known_indels_vcf,
        reference = reference,
        runTimeSettings = runTimeSettings
    }

    call SNPGenotyping.SNPGenotyping {
      input:
        sample_id = sample_id,
        output_basename = output_basename,
        input_bam = ImportShortReadAlignment.output_bam,
        input_bam_index = ImportShortReadAlignment.output_bam_index,
        alleles_vcf = alleles_vcf,
        alleles_vcf_index = alleles_vcf_index,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  output {
    Array[File] indel_realigner_interval_list_file = ImportShortReadAlignment.indel_realigner_interval_list_file
    Array[File] output_bam = ImportShortReadAlignment.output_bam
    Array[File] output_bam_index = ImportShortReadAlignment.output_bam_index
    Array[File] validate_samfile_report_file = ImportShortReadAlignment.validate_samfile_report_file
    Array[File] samtools_stats_report_file = ImportShortReadAlignment.samtools_stats_report_file
    Array[File] samtools_idxstats_report_file = ImportShortReadAlignment.samtools_idxstats_report_file
    Array[File] samtools_flagstat_report_file = ImportShortReadAlignment.samtools_flagstat_report_file
    Array[File] callable_loci_summary_file = ImportShortReadAlignment.callable_loci_summary_file
    Array[File] output_vcf = SNPGenotyping.output_vcf
    Array[File] output_vcf_index = SNPGenotyping.output_vcf_index
    Array[File] zarr_output = SNPGenotyping.zarr_output
  }
}
