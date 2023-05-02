version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Target Regions portion of the vector Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.

import "../../../structs/gcp/RunTimeSettings.wdl"
import "../../../tasks/gcp/Tasks.wdl" as Tasks

workflow TargetRegions {
  String pipeline_version = "1.0.0"

  input {
    String project_id
    String sample_name
    File input_bam
    File input_bam_index
    File sample_manifest
    File gene_coordinates_file
    File metadata
    File species_id_file

    File CNV_HMM_output # zip of the coveragefolder
    File HMM_coverage_variance_file
    File plotting_functions_file

    RunTimeSettings runTimeSettings
    String runtime_zones
  }

  String output_basename = project_id + "_" + sample_name

  # Step 1: Extract Diagnostic Reads
  call TargetRegionsTasks.ExtractDiagnosticReads as ExtractDiagnosticReads {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      project_id = project_id,
      output_basename = output_basename,
      runTimeSettings = runTimeSettings,
      runtime_zones = runtime_zones
  }

  # Step 2: Target Regions CNV calling
  call TargetRegionsTasks.TargetRegionsCNVCalling as CNVCalling {
    input:
      sample_manifest = sample_manifest
      gene_coordinates_file
      sample_metadata
      species_id_file
      coverage_variance_file
      coverage_tar
      diagnostic_reads_tar
      plotting_functions_file

      HMM = CNV_HMM_output,
      discordant_reads = ExtractDiagnosticReads.discordant_reads_output,
      breakpoint_reads = ExtractDiagnosticReads.breakpoint_reads_output,
      output_basename = output_basename,
      runTimeSettings = runTimeSettings,
      runtime_zones = runtime_zones
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File output_file = CNVCalling.output_file
  }
}

task ExtractDiagnosticReads {
  input {
    File input_bam
    File input_bam_index
    String project_id
    String output_basename

    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
    Int preemptible_tries = runTimeSettings.preemptible_tries
    String runtime_zones = "us-central1-b"
    Float mem_gb = 3.75
    Int disk_gb = 50
  }

  command <<<
    # Get the discordant reads
    # Runs SSFA.py for every chromosome: This script goes through an alignment file and records the positions of reads within a specified region whose mates map to a different chromosome or discordantly on the same chromosome
    SSFA_script=$scriptsfolder/SSFA.py
    SSFAfolder=$outputfolder/diagnostic_reads/SSFA
    python $SSFA_script $bamfile 2R 3425000:3650000 ${SSFAfolder}/2R/Ace1_region/${samplename}_Ace1_SSFA_output.csv 10 > ${SSFAfolder}/2R/Ace1_region/SSFAlogs/${samplename}_Ace1_SSFA_output.log 2>&1
    python $SSFA_script $bamfile 2R 28460000:28570000 ${SSFAfolder}/2R/Cyp6_region/${samplename}_CYP6_SSFA_output.csv 10 > ${SSFAfolder}/2R/Cyp6_region/SSFAlogs/${samplename}_CYP6_SSFA_output.log 2>&1
    python $SSFA_script $bamfile 3R 6900000:7000000 ${SSFAfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_SSFA_output.csv 10 > ${SSFAfolder}/3R/Cyp6zm_region/SSFAlogs/${samplename}_CYP6ZM_SSFA_output.log 2>&1
    python $SSFA_script $bamfile 3R 28570000:28620000 ${SSFAfolder}/3R/Gste_region/${samplename}_GST_SSFA_output.csv 10 > ${SSFAfolder}/3R/Gste_region/SSFAlogs/${samplename}_GST_SSFA_output.log 2>&1
    python $SSFA_script $bamfile X 15220000:15255000 ${SSFAfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_SSFA_output.csv 10 > ${SSFAfolder}/X/Cyp9k1_region/SSFAlogs/${samplename}_CYP9K1_SSFA_output.log 2>&1

    # Get the soft clipped reads
    # Runs breakpoint_detector.py for every chrom: This script goes through an alignment file and records the positions at which soft_clipping is detected in the aligned reads
    breakpoints_script=$scriptsfolder/breakpoint_detector.py
    breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
    python $breakpoints_script $bamfile 2R 3425000:3650000 ${breakpointsfolder}/2R/Ace1_region/${samplename}_Ace1_breakpoints_output 10 > ${breakpointsfolder}/2R/Ace1_region/breakpointlogs/${samplename}_Ace1_breakpoints_output.log 2>&1
    python $breakpoints_script $bamfile 2R 28460000:28570000 ${breakpointsfolder}/2R/Cyp6_region/${samplename}_CYP6_breakpoints_output 10 > ${breakpointsfolder}/2R/Cyp6_region/breakpointlogs/${samplename}_CYP6_breakpoints_output.log 2>&1
    python $breakpoints_script $bamfile 3R 6900000:7000000 ${breakpointsfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_breakpoints_output 10 > ${breakpointsfolder}/3R/Cyp6zm_region/breakpointlogs/${samplename}_CYP6ZM_breakpoints_output.log 2>&1
    python $breakpoints_script $bamfile 3R 28570000:28620000 ${breakpointsfolder}/3R/Gste_region/${samplename}_GST_breakpoints_output 10 > ${breakpointsfolder}/3R/Gste_region/breakpointlogs/${samplename}_GST_breakpoints_output.log 2>&1
    python $breakpoints_script $bamfile X 15220000:15255000 ${breakpointsfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_breakpoints_output 10 > ${breakpointsfolder}/X/Cyp9k1_region/breakpointlogs/${samplename}_CYP9K1_breakpoints_output.log 2>&1

  >>>
  runtime {
    docker: docker
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: mem_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    zones: runtime_zones
  }
  output {
    File diagnostic_reads_tar = # zip of the directory structure for the extract diagnostic reads step
  }
}

task TargetRegionsCNVCalling {
  input {
    File sample_manifest              # manifest: pipeline input
    File gene_coordinates_file        # gene_coordinates_file: pipeline input
    File sample_metadata              # metadata: pipeline input
    File species_id_file              # species_id_file: pipeine input
    File coverage_variance_file       # coverage_variance_file: output from CoverageSummary step in HMM pipeline
    File coverage_tar                 # coveragefolder: output from CoverageSummary step in HMM pipeline
    File diagnostic_reads_tar         # diagnostic_reads_folder: output of ExtractDiagnosticReads
    File plotting_functions_file      # plotting_functions_file

    String project_id
    String output_basename

    String docker = "us.gcr.io/broad-gotc-prod/r:3.6.1"
    Int num_cpu = 8
    RunTimeSettings runTimeSettings
    Int preemptible_tries = runTimeSettings.preemptible_tries
    String runtime_zones = "us-central1-b"
    Float mem_gb = 3.75
    Int disk_gb = 50
  }


  command <<<
    # need to unzip the tarred folders before passing them here
    R-3.6.1 --slave -f $scriptsfolder/target_regions_analysis.r --args ~{sample_manifest}# $manifest \
      ~{gene_coordinates_file} \      # $gene_coordinates_file
      ~{sample_metadata} \            # $metadata
      ~{species_id_file} \            # $species_id_file
      ~{coverage_variance_file} \     # $coverage_variance_file
      ~{coverage_tar} \               # $coveragefolder
      ~{diagnostic_reads_tar} \       # $diagnostic_reads_folder
      ~{plotting_functions_file} \    # $plotting_functions_file
      ~{num_cpu} \                    # $ncores
      > target_regions_analysis/target_regions_analysis.log 2>&1
  >>>
  runtime {
    docker: docker
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: mem_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    zones: runtime_zones
  }
  output {
    File output_file = local_file
    File output_index_file = "~{local_file}.tbi"
  }
}
