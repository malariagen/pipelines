version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2021
##
## This WDL pipeline implements the Read-Backed component of the Mospquito Phasing Pipeline as described in
## https://github.com/malariagen/pipelines/blob/master/docs/specs/phasing-vector.md
##


workflow HMM {
  String pipeline_version = "1.0.0"

  input {
    File input_bam
    String sample_name
    String output_dir
  }

  call WindowedCoverage {
    input:
      input_bam = input_bam,
      sample_name = sample_name,
      output_dir = output_dir
  }

  output {
    String out = ""
  }
}

task WindowedCoverage {
  input {
    File input_bam
    String sample_name
    String output_dir

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    Int ram = "8 GiB"
    Int cpu = 16
    # TODO: Make disk space dynamic based on input size
    Int disk = 70
    Int preemptible = 3
  }
  command <<<
    echo "Processing file: " 
    basename ~{input_bam}
    echo "Current directory: " 
    pwd
    cd /cnv/scripts/
    ls -lht
    bash get_windowed_coverage_and_diagnostic_reads.sh ~{input_bam} ~{sample_name} ~{output_dir}
  >>>
  runtime {
    docker: docker
    memory: "${ram} GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  output {
    String message = read_string(stdout())
  }
}


