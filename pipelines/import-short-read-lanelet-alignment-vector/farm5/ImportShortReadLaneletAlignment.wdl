version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
## Imports multiple lanelet bams/crams for the same sample from IRODS.
## BWA aligns the lanelet bams/crams.  Processes each lanelet in parallel.
## Lanelet refers a lane of a sample that has been multiplexed into multiple
## sequencing lanes.
## Expects that the per_sample_manifest_file contains only the lanelets for a single sample
## and each lanelet has its own row.
## Columns:  sample, run_ena, irods_path

import "../../../structs/farm5/RunTimeSettings.wdl"
import "../../../structs/ReferenceSequence.wdl"
import "../../../tasks/farm5/ImportTasks.wdl" as ImportTask
import "../../../tasks/farm5/Alignment.wdl" as LaneletAlignmentTask
import "../../../tasks/farm5/ShortReadAlignmentTasks.wdl" as AlignTasks



workflow ImportShortReadLaneletAlignment {
  String pipeline_version = "1.0.0"
  
  String LANELET_INFO_COLNAME_IRODS_PATH = "irods_path"

  input {
    String sample_id
    File per_sample_manifest_file

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Array[Object] lanelet_infos = read_objects(per_sample_manifest_file)
  scatter (idx in range(length(lanelet_infos))) {

    String irods_path = lanelet_infos[idx][LANELET_INFO_COLNAME_IRODS_PATH]

    call ImportTask.ImportIRODS as ImportIRODS {
      input:
        irods_path = irods_path,
        sample_id = sample_id,
        runTimeSettings = runTimeSettings
    }

    Boolean is_cram = sub(basename(irods_path), ".*\\.", "") == "cram"
    String cram_path = if is_cram then ImportIRODS.output_file else ""
    String bam_path = if is_cram then "" else ImportIRODS.output_file

    String irods_basename = if is_cram then basename(irods_path, ".cram")  else basename(irods_path, ".bam")
    String read_group_filename = sample_id + "." + irods_basename + ".read_group.txt"

    call AlignTasks.ExtractReadGroup as ExtractReadGroup {
      input:
        input_file = ImportIRODS.output_file,
        sample_id = sample_id,
        read_group_filename = read_group_filename,
        runTimeSettings = runTimeSettings
    }

    String read_group = read_string(ExtractReadGroup.read_group_file)
    String align_basename = sample_id + "." + irods_basename + ".align"

    call LaneletAlignmentTask.Alignment as LaneletAlignment {
      input:
        sample_id = sample_id,
        read_group = read_group,
        input_cram = cram_path,
        input_bam = bam_path,
        input_fastq1 = "",
        input_fastq2 = "",
        output_file_basename = align_basename,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  output {
    Array[File] output_bam = LaneletAlignment.output_bam
  }

  meta {
    allowNestedInputs: true
  }
}
