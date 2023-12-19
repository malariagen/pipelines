version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements CNV Coverage Calls portion of the vector Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.

workflow CNVCoverageCalls {
  String pipeline_version = "1.0.0"

  ## Sub-pipeline: CNV Coverage Calls

  input {
    String chromosome               # chrom: Runs separately for each species and chromosome (2L, 2R, 3L, 3R, X)
    File sample_species_manifest    # manifest: species specific manifest file - NOTE different from other manifests
    File gene_coordinates_file      # gene_coordinates_file: pipeline input
    File detox_genes_file           # detox_genes_file: pipeline input
    File consolidated_coverage_dir_tar # workingfolder: Specifies the folder containing the HMM output files - will be a tarred output here
    File sample_metadata            # metadata: pipeline input
    String species                  # species: pipeline input
    Int num_samples                 # number of samples - used to calucluate the memory and disk requested for coverage calls step
    Int preemptible_tries
    String runtime_zones
  }

  # Step 1: CNV Coverage Calls
  call CoverageCalls {
    input:
      chromosome = chromosome,
      sample_species_manifest = sample_species_manifest,
      gene_coordinates_file = gene_coordinates_file,
      detox_genes_file = detox_genes_file,
      consolidated_coverage_dir_tar = consolidated_coverage_dir_tar,
      sample_metadata = sample_metadata,
      species = species,
      num_samples = num_samples,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File cnv_coverage_table = CoverageCalls.full_coverage_CNV_table
    File cnv_coverage_raw_table = CoverageCalls.full_raw_CNV_table
    File cnv_coverage_Rdata = CoverageCalls.coverage_CNV_Rdata
  }
}

task CoverageCalls {
  input {
    String chromosome               # chrom: Runs separately for each species and chromosome (2L, 2R, 3L, 3R, X)
    File sample_species_manifest    # manifest: species specific manifest file - NOTE different from other manifests
    File gene_coordinates_file      # gene_coordinates_file: pipeline input
    File detox_genes_file           # detox_genes_file: pipeline input
    File consolidated_coverage_dir_tar # workingfolder: Specifies the folder containing the HMM output files - will be a tarred output here
    File sample_metadata            # metadata: pipeline input
    String species                  # species: pipeline input
    Int num_samples                 # number of samples - used to calucluate the memory and disk requested for coverage calls step

    # Runtime Settings
    Int preemptible_tries
    String runtime_zones
    String docker = "us.gcr.io/broad-gotc-prod/cnv/r:1.0.0-1692386236"
    Int num_cpu = 1
    Int memory_gb = if (floor(num_samples * 0.2) + 3.75) < 120 then (floor(num_samples * 0.2) + 3.75) else 120
    Int disk_size_gb = if (ceil(num_samples * 0.75) + 10) < 1200 then (ceil(num_samples * 0.75) + 10) else 1200
  }
  # ncores: number of CPUs - has not been optimized
  String output_dir = species + "_CNV"

  # These will be named by the script and cannot be changed here
  String full_coverage_CNV_table = output_dir + '/full_coverage_CNV_table_' + chromosome + '.csv'
  String full_raw_CNV_table = output_dir + '/full_raw_CNV_table_'  + chromosome + '.csv'
  String Rdata = output_dir + '/CNV_analysis_' + chromosome +'.Rdata'
  # Since I am not outputting the full directory with the species included in the folder structure, I am renaming the output file to include the species name
  String species_full_coverage_CNV_table = output_dir + '/full_coverage_CNV_table_' + species + '_' + chromosome + '.csv'
  String species_full_raw_CNV_table = output_dir + '/full_raw_CNV_table_' + species + '_' + chromosome + '.csv'
  String species_Rdata = output_dir + '/CNV_analysis_' + species + '_' + chromosome +'.Rdata'

  # once we untarr the HMM dir, the path to the untarred dir will be local to where the command runs
  String coverage_dir = basename(consolidated_coverage_dir_tar, ".tar.gz")
  String HMM_working_dir = coverage_dir + '/' + chromosome + '/HMM_output'



  command <<<
    # unzip the coverage tarball and get the directory name
    tar -zxvf ~{consolidated_coverage_dir_tar}

    # set up the output dir
    mkdir ~{HMM_working_dir}/~{output_dir}

    # The following combines the coverage_variance files into a signle file
    # Get the name of the first coverage variance file - doesn't matter which one
    FIRST_COVERAGE_VARIANCE_FILE=$(ls ~{coverage_dir}/coverage_variance_masked_* | head -n 1)

    # Add the header to the new combined file (the header is in every file, but we only want it once in the combined file)
    HEADER=$(head -n 1 $FIRST_COVERAGE_VARIANCE_FILE)
    echo "sample_id$HEADER" >  single_coverage_variance_masked.csv

    # Skip the header and add the remaining contents of each coverage variance file to the combined file
    for f in $(ls ~{coverage_dir}/coverage_variance_masked_*); do grep -v "$HEADER" $f >> single_coverage_variance_masked.csv; done

    echo "single_coverage_variance_masked.csv"
    cat single_coverage_variance_masked.csv

    echo "Starting R script"
    /opt/R/3.6.1/bin/R --slave -f /usr/local/Rscripts/CNV_analysis.r --args ~{chromosome} \
      ~{sample_species_manifest} \
      single_coverage_variance_masked.csv \
      ~{gene_coordinates_file} \
      ~{detox_genes_file} \
      ~{HMM_working_dir} \
      ~{num_cpu} \
      ~{output_dir}\
      ~{sample_metadata}
    echo "R script complete"

    mv ~{HMM_working_dir}/~{full_coverage_CNV_table} ~{species_full_coverage_CNV_table}
    mv ~{HMM_working_dir}/~{full_raw_CNV_table} ~{species_full_raw_CNV_table}
    mv ~{HMM_working_dir}/~{Rdata} ~{species_Rdata}
  >>>
  runtime {
    docker: docker
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_size_gb + " HDD"
    zones: runtime_zones
  }
  output {
    File full_coverage_CNV_table = species_full_coverage_CNV_table
    File full_raw_CNV_table = species_full_raw_CNV_table
    File coverage_CNV_Rdata = species_Rdata
  }
}
