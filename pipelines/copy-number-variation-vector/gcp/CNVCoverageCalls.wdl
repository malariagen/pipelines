version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements CNV Coverage Calls portion of the vector Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.

import "../../../tasks/gcp/Tasks.wdl" as Tasks
import "../../../tasks/gcp/StatisticalPhasingTasks.wdl" as StatisticalPhasingTasks

workflow CNVCoverageCalls {
  String pipeline_version = "1.0.0"

  ## Sub-pipeline: CNV Coverage Calls

  input {
    String project_id
    String sample_name
    File input_bam
    File input_bam_index

    File CNV_HMM_output

    String chromosome
    File sample_species_manifest
    File coverage_variance_file
    File gene_coordinates_file
    File detox_genes_file
    File HMM_working_dir_tar
    File sample_metadata

    Int preemptible_tries
    String runtime_zones
  }

  String output_basename = project_id + "_" + sample_name

  # Step 1: CNV Coverage Calls
  call CNVCoverageTasks.CNVCoverageCalls as CoverageCalls {
    input:
    chromosome = chromosome,
    sample_species_manifest = sample_species_manifest,
    coverage_variance_file = coverage_variance_file,
    gene_coordinates_file = gene_coordinates_file,
    detox_genes_file = detox_genes_file,
    HMM_working_dir_tar = HMM_working_dir_tar,
    sample_metadata = sample_metadata,
    preemptible_tries = preemptible_tries,
    runtime_zones = runtime_zones,
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File cnv_coverage_csv = CoverageCalls.cnv_coverage_csv
    File cnv_coverage_Rdata = CoverageCalls.cnv_coverage_Rdata
  }
}

task CNVCoverageCalls {
  input {
    String chromosome               # chrom: Runs separately for each species and chromosome (2L, 2R, 3L, 3R, X)
    File sample_species_manifest    # manifest: species specific manifest file - NOTE different from other manifests
    File coverage_variance_file     # coverage_variance_file: output from CoverageSummary step in HMM pipeline
    File gene_coordinates_file      # gene_coordinates_file: pipeline input
    File detox_genes_file           # detox_genes_file: pipeline input
    File HMM_working_dir_tar        # workingfolder: Specifies the folder containing the HMM output files - will be a tarred output here
    File sample_metadata            # metadata: pipeline input
    String project_id
    String species

    # Runtime Settings
    Int preemptible_tries
    String runtime_zones
    String docker = "us.gcr.io/broad-gotc-prod/r:3.6.1"
    Int num_cpu = 1
    Float mem_gb = 3.75
    Int disk_gb = 50
  }
  # ncores: number of CPUs - has not been optimized
  # outputfolder: Make sure that the output name contains the name of the species (since this will be run once per species)
    String output_dir = project_id + "_" + species
    String output_name = species + "_CNV"
    # coveragefolder example: /lustre/scratch118/malaria/team112/personal/el10/v3.7_1246-VO-TZ-KABULA-VMF00185/coverage
    # here I have just used "coverage"
    String full_coverage_CNV_table = output_dir + '/full_coverage_CNV_table_' + chromosome + '.csv'
    String full_raw_CNV_table = output_dir + '/full_raw_CNV_table_'+ chromosome + '.csv'
    String Rdata = output_dir + '/CNV_analysis_' + chromosome +'.Rdata'

  command <<<
    # set up directories as needed

    R-3.6.1 --slave -f $scriptsfolder/CNV_analysis.r --args ~{chromosome} \
      ~{sample_species_manifest} \
      ~{coverage_variance_file} \
      ~{gene_coordinates_file} \
      ~{detox_genes_file} \
      ~{HMM_working_dir_tar} \
      ~{num_cpu} \
      ~{output_dir}\
      ~{sample_metadata} \
      > coverage/~{chromosome}/CNV_analysis_logs/CNV_analysis_${output_name}.log 2>&1
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
    File full_coverage_CNV_table = full_coverage_CNV_table
    File full_raw_CNV_table = full_raw_CNV_table
    File coverage_CNV_Rdata = Rdata
  }
}
