version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.

import "HMM.wdl" as HMM
import "TargetRegions.wdl" as TargetRegions
import "CNVCoverageCalls.wdl" as CNVCoverageCalls
import "CNVTasks.wdl" as CNVTasks

workflow CNV {
  meta {
    description: "This is a pipeline for calling Copy Number Variants (CNVs) for a cohort of multiple samples. This inclludes an HMM step followed by coverage calls and target regions pipelines to improve acuracy."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    # windowed coverage inputs
    Array[File] input_bams
    Array[File] input_bais
    Array[String] sample_ids
    String output_dir="coverage"
    Int interval = 300
    Int window_size = 300
    Int min_qual = 10
    # coverage summary stats and coverage HMM inputs
    Float accessibility_threshold = 0.9
    Float mapq_threshold = 0.5
    File accessibility_mask_file
    File mapq_file
    File sample_manifest
    File gc_content_file
    String sample_group_id
    # coverage HMM inputs
    String species
    # target regions inputs
    File gene_coordinates_file
    File sample_metadata
    File species_id_file
    File plotting_functions_file
    # coverage calls inputs
    Array[String] chromosomes = ["2L", "2R", "3L", "3R", "X"]
    File detox_genes_file
    # runtime inputs
    Int preemptible_tries
    String runtime_zones = "us-central1-b"
  }

  # This is a wdl hack to create a pseudo None
  if (false) {
    File? none = "None"
  }

  # Scatter over samples
  scatter(idx in range(length(input_bams))) {
    # Run windowed coverage on each sample (for each chromosome)
    call HMM.HMM as HMM {
      input:
        input_bam = input_bams[idx],
        input_bai = input_bais[idx],
        sample_id = sample_ids[idx],
        output_dir = output_dir,
        interval = interval,
        window_size = window_size,
        min_qual = min_qual,
        accessibility_threshold = accessibility_threshold,
        mapq_threshold = mapq_threshold,
        accessibility_mask_file = accessibility_mask_file,
        mapq_file = mapq_file,
        sample_manifest = sample_manifest,
        gc_content_file = gc_content_file,
        sample_group_id = sample_group_id,
        species = species,
        runtime_zones = runtime_zones
    }

    call TargetRegions.TargetRegions as TargetRegions {
      input:
        sample_id = sample_ids[idx],
        input_bam = input_bams[idx],
        input_bam_index = input_bais[idx],
        sample_manifest = sample_manifest,
        gene_coordinates_file = gene_coordinates_file,
        sample_metadata = sample_metadata,
        species_id_file = species_id_file,
        CNV_HMM_output = HMM.output_gz, # zip of the coveragefolder
        HMM_coverage_variance_file = HMM.coverage_variance,
        plotting_functions_file = plotting_functions_file,
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones
    }
  }

  call CNVTasks.ConsolidateHMMOutput as CHMM {
    input:
      hmm_tarballs = HMM.output_gz,
      output_dir = output_dir,
      runtime_zones = runtime_zones
  }

  call CNVTasks.CreateSampleManifest as CreateSampleManifest {
    input:
      sample_ids = sample_ids,
      runtime_zones = runtime_zones
  }

  # Runs separately for each chromosome (2L, 2R, 3L, 3R, X)
  scatter(idx in range(length(chromosomes))) {
    call CNVCoverageCalls.CNVCoverageCalls as CNVCoverageCalls {
      input:
        chromosome = chromosomes[idx],
        sample_species_manifest = CreateSampleManifest.sample_manifest,
        gene_coordinates_file = gene_coordinates_file,
        detox_genes_file = detox_genes_file,
        consolidated_coverage_dir_tar = CHMM.consolidated_gz,
        sample_metadata = sample_metadata,
        species = species,
        num_samples = length(sample_ids),
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones
    }
  }

  output {
    File hmm_tar = CHMM.consolidated_gz
    Array[File] targeted_regions_focal_region_cnv_table = TargetRegions.focal_region_CNV_table
    Array[File] targeted_regions_hmm_gene_copy_number = TargetRegions.HMM_gene_copy_number
    Array[File] cnv_coverage_calls_tables = CNVCoverageCalls.cnv_coverage_table
    Array[File] cnv_coverage_calls_raw_tables = CNVCoverageCalls.cnv_coverage_raw_table
    Array[File] cnv_coverage_calls_Rdata = CNVCoverageCalls.cnv_coverage_Rdata
  }
}