version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.

import "HMM.wdl" as HMM
import "TargetRegions.wdl" as TargetRegions

workflow CNV {
  meta {
    description: "This is a pipeline for calling Copy Number Variants (CNVs) for a cohort of multiple samples. This inclludes an HMM step followed by coverage calls and target regions pipelines to improve acuracy."
    allowNestedInputs: true
  }

  String pipeline_version = "1.0.0"

  input {
    String project_id
    # windowed coverage inputs
    Array[File] input_bams
    Array[File] input_bais
    Array[String] sample_names
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
    #coverage HMM inputs
    String species

    File gene_coordinates_file
    File sample_metadata
    File species_id_file
    File plotting_functions_file

    Int preemptible_tries
    String runtime_zones
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
        sample_name = sample_names[idx],
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
        species = species
    }

    call TargetRegions.TargetRegions as TargetRegions {
      input:
        sample_name = sample_names[idx],
        input_bam = input_bams[idx],
        input_bam_index = input_bais[idx],
        sample_manifest = sample_manifest,
        gene_coordinates_file = gene_coordinates_file,
        sample_metadata = sample_metadata,
        species_id_file = species_id_file,
        CNV_HMM_output = HMM.output_gz ,# zip of the coveragefolder
        HMM_coverage_variance_file = HMM.coverage_variance ,
        plotting_functions_file = plotting_functions_file,
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones
    }
  }

  call ConsolidateHMMOutput as CHMM {
    input:
      hmm_tarballs = HMM.output_gz,
      output_dir = output_dir
  }

  output {
    File hmm_tar = CHMM.consolidated_gz
    Array[File] targeted_regions_focal_region_cnv_table = TargetRegions.focal_region_CNV_table
    Array[File] targeted_regions_hmm_gene_copy_number = TargetRegions.HMM_gene_copy_number
  }
}

task ConsolidateHMMOutput {
  meta {
    description: "This task takes the output from the HMM subpipeline and combines it into a single tarball."
  }
  parameter_meta {
    hmm_tarballs: "The output files from the HMM sub-pipeline. This is an array of tarballs, one for each sample."
    output_dir: "The output directory for the tarball."
  }
  input {
    Array[File] hmm_tarballs
    String output_dir
    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    # TODO: Make disk space dynamic based on input size
    Int disk = 70
    Int preemptible = 3
  }
  command <<<
    set -x
    echo "Current directory: "
    pwd
    #For each file in hmm_tarballs, extract the tarball and move the contents to the output directory
    for tarball in ~{sep=' ' hmm_tarballs} ; do
      echo "Extracting tarball: " $tarball
      tar --backup=numbered -zxvf $tarball
      ls -lht
    done
    ls -lht
    tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
  >>>
  runtime {
    docker: docker
    memory: ram
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File consolidated_gz = "~{output_dir}.tar.gz"
  }
}