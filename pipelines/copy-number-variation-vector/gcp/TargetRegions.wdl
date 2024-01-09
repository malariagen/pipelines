version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2023
##
## This WDL pipeline implements Target Regions portion of the vector Copy Number Variation Pipeline as described in:
## https://github.com/malariagen/pipelines/blob/add_cnv_vector_spec/docs/specs/cnv-vector.md
##.


workflow TargetRegions {
  String pipeline_version = "1.0.0"

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File gene_coordinates_file
    File sample_metadata
    File species_id_file

    File CNV_HMM_output # zip of the coveragefolder
    File HMM_coverage_variance_file
    File plotting_functions_file

    Int preemptible_tries
    String runtime_zones
  }

  # Step 1: Extract Diagnostic Reads
  call ExtractDiagnosticReads as ExtractDiagnosticReads {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      sample_id = sample_id,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  # Step 2: Target Regions CNV calling
  call TargetRegionsCNVCalling as CNVCalling {
    input:
      sample_id = sample_id,
      gene_coordinates_file = gene_coordinates_file,
      sample_metadata = sample_metadata,
      species_id_file = species_id_file,
      coverage_variance_file = HMM_coverage_variance_file,
      coverage_tar = CNV_HMM_output,
      diagnostic_reads_tar = ExtractDiagnosticReads.diagnostic_reads_tar,
      plotting_functions_file = plotting_functions_file,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  meta {
    allowNestedInputs: true
  }

  output {
    File diagnostic_reads_tar = ExtractDiagnosticReads.diagnostic_reads_tar
    File focal_region_CNV_table = CNVCalling.focal_region_CNV_table
    File HMM_gene_copy_number = CNVCalling.HMM_gene_copy_number
    #File target_regions_Rdata = CNVCalling.target_regions_Rdata
  }
}

task ExtractDiagnosticReads {
  input {
    File input_bam
    File input_bam_index
    String sample_id

    Int preemptible_tries
    String runtime_zones
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    Int num_cpu = 1
    Float mem_gb = 3.75
    Int disk_gb = 50
  }


  command <<<
    set -euo pipefail
    # Get the discordant reads
    # Runs SSFA.py for every chromosome: This script goes through an alignment file and records the positions of reads within a specified region whose mates map to a different chromosome or discordantly on the same chromosome

    # create directories for diagnostic reads output files and logs
    mkdir -p diagnostic_reads/SSFA/2R/Ace1_region/SSFAlogs
    mkdir -p diagnostic_reads/SSFA/2R/Cyp6_region/SSFAlogs
    mkdir -p diagnostic_reads/SSFA/3R/Cyp6zm_region/SSFAlogs
    mkdir -p diagnostic_reads/SSFA/3R/Gste_region/SSFAlogs
    mkdir -p diagnostic_reads/SSFA/X/Cyp9k1_region/SSFAlogs

    # create directories for diagnostic reads output files and logs
    mkdir -p diagnostic_reads/breakpoints/2R/Ace1_region/breakpointlogs
    mkdir -p diagnostic_reads/breakpoints/2R/Cyp6_region/breakpointlogs
    mkdir -p diagnostic_reads/breakpoints/3R/Cyp6zm_region/breakpointlogs
    mkdir -p diagnostic_reads/breakpoints/3R/Gste_region/breakpointlogs
    mkdir -p diagnostic_reads/breakpoints/X/Cyp9k1_region/breakpointlogs

    #activate the conda environment
    source activate cnv37

    python /cnv/scripts/SSFA.py ~{input_bam} 2R 3425000:3650000 diagnostic_reads/SSFA/2R/Ace1_region/~{sample_id}_Ace1_SSFA_output.csv 10 > diagnostic_reads/SSFA/2R/Ace1_region/SSFAlogs/~{sample_id}_Ace1_SSFA_output.log 2>&1
    python /cnv/scripts/SSFA.py ~{input_bam} 2R 28460000:28570000 diagnostic_reads/SSFA/2R/Cyp6_region/~{sample_id}_CYP6_SSFA_output.csv 10 > diagnostic_reads/SSFA/2R/Cyp6_region/SSFAlogs/~{sample_id}_CYP6_SSFA_output.log 2>&1
    python /cnv/scripts/SSFA.py ~{input_bam} 3R 6900000:7000000 diagnostic_reads/SSFA/3R/Cyp6zm_region/~{sample_id}_CYP6ZM_SSFA_output.csv 10 > diagnostic_reads/SSFA/3R/Cyp6zm_region/SSFAlogs/~{sample_id}_CYP6ZM_SSFA_output.log 2>&1
    python /cnv/scripts/SSFA.py ~{input_bam} 3R 28570000:28620000 diagnostic_reads/SSFA/3R/Gste_region/~{sample_id}_GST_SSFA_output.csv 10 > diagnostic_reads/SSFA/3R/Gste_region/SSFAlogs/~{sample_id}_GST_SSFA_output.log 2>&1
    python /cnv/scripts/SSFA.py ~{input_bam} X 15220000:15255000 diagnostic_reads/SSFA/X/Cyp9k1_region/~{sample_id}_CYP9K1_SSFA_output.csv 10 > diagnostic_reads/SSFA/X/Cyp9k1_region/SSFAlogs/~{sample_id}_CYP9K1_SSFA_output.log 2>&1

    # Get the soft clipped reads
    # Runs breakpoint_detector.py for every chrom: This script goes through an alignment file and records the positions at which soft_clipping is detected in the aligned reads

    python /cnv/scripts/breakpoint_detector.py ~{input_bam} 2R 3425000:3650000 diagnostic_reads/breakpoints/2R/Ace1_region/~{sample_id}_Ace1_breakpoints_output 10 > diagnostic_reads/breakpoints/2R/Ace1_region/breakpointlogs/~{sample_id}_Ace1_breakpoints_output.log 2>&1
    python /cnv/scripts/breakpoint_detector.py ~{input_bam} 2R 28460000:28570000 diagnostic_reads/breakpoints/2R/Cyp6_region/~{sample_id}_CYP6_breakpoints_output 10 > diagnostic_reads/breakpoints/2R/Cyp6_region/breakpointlogs/~{sample_id}_CYP6_breakpoints_output.log 2>&1
    python /cnv/scripts/breakpoint_detector.py ~{input_bam} 3R 6900000:7000000 diagnostic_reads/breakpoints/3R/Cyp6zm_region/~{sample_id}_CYP6ZM_breakpoints_output 10 > diagnostic_reads/breakpoints/3R/Cyp6zm_region/breakpointlogs/~{sample_id}_CYP6ZM_breakpoints_output.log 2>&1
    python /cnv/scripts/breakpoint_detector.py ~{input_bam} 3R 28570000:28620000 diagnostic_reads/breakpoints/3R/Gste_region/~{sample_id}_GST_breakpoints_output 10 > diagnostic_reads/breakpoints/3R/Gste_region/breakpointlogs/~{sample_id}_GST_breakpoints_output.log 2>&1
    python /cnv/scripts/breakpoint_detector.py ~{input_bam} X 15220000:15255000 diagnostic_reads/breakpoints/X/Cyp9k1_region/~{sample_id}_CYP9K1_breakpoints_output 10 > diagnostic_reads/breakpoints/X/Cyp9k1_region/breakpointlogs/~{sample_id}_CYP9K1_breakpoints_output.log 2>&1

    ls -lht
    tar -zcvf diagnostic_reads.tar.gz diagnostic_reads
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
    File diagnostic_reads_tar = "diagnostic_reads.tar.gz"
  }
}

task TargetRegionsCNVCalling {
  input {
    String sample_id
    File gene_coordinates_file        # gene_coordinates_file: pipeline input
    File sample_metadata              # metadata: pipeline input
    File species_id_file              # species_id_file: pipeine input
    File coverage_variance_file       # coverage_variance_file: output from CoverageSummary step in HMM pipeline
    File coverage_tar                 # coveragefolder: output from CoverageSummary step in HMM pipeline
    File diagnostic_reads_tar         # diagnostic_reads_folder: output of ExtractDiagnosticReads
    File plotting_functions_file      # plotting_functions_file

    String docker = "us.gcr.io/broad-gotc-prod/cnv/r:1.0.0-1692386236"
    Int num_cpu = 8
    Int preemptible_tries
    String runtime_zones = "us-central1-b"
    Float mem_gb = 3.75
    Int disk_gb = 50
  }

  String coverage_dir = basename(coverage_tar, ".tar.gz")
  String diagnostic_reads_dir = basename(diagnostic_reads_tar, ".tar.gz")

  command <<<
    set -euo pipefail

    # get the manifest, metadata, and species_id files for just this sample
    echo "~{sample_id}" > ~{sample_id}_manifest.tsv

    head -n 1 ~{sample_metadata} >> ~{sample_id}_metadata.tsv
    grep "~{sample_id}" ~{sample_metadata} >> ~{sample_id}_metadata.tsv

    head -n 1 ~{species_id_file} >> ~{sample_id}_species_id.tsv
    grep "~{sample_id}" ~{species_id_file} >> ~{sample_id}_species_id.tsv

    # unzip the coverage tarball and get the directory name
    tar -zxvf ~{coverage_tar}

    # unzip the diagnostic reads tarball and get the directory name
    tar -zxvf ~{diagnostic_reads_tar}

    #create output directory
    mkdir target_regions_analysis

    echo "Starting R script"
    /opt/R/3.6.1/bin/R --slave -f /usr/local/Rscripts/target_regions_analysis.r --args ~{sample_id}_manifest.tsv \
      ~{gene_coordinates_file} \
      ~{sample_id}_metadata.tsv \
      ~{sample_id}_species_id.tsv \
      ~{coverage_variance_file} \
      ~{coverage_dir} \
      ~{diagnostic_reads_dir} \
      ~{plotting_functions_file} \
      ~{num_cpu}
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
    File focal_region_CNV_table = "target_regions_analysis/focal_region_CNV_table.csv"
    File HMM_gene_copy_number = "target_regions_analysis/HMM_gene_copy_number.csv"
    # File target_regions_Rdata = "target_regions_analysis/target_regions_analysis.Rdata"
  }
}
