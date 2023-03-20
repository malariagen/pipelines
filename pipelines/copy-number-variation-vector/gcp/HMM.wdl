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
  }

  call WindowedCoverage {
    input:
      input_bam = input_bam
  }

  output {
    String out = ""
  }
}

task WindowedCoverage {
  input {
    File input_bam
  }
  command {
    echo "Hello, World!"
    basename ~{input_bam}
    pwd
    cd scripts/
    bash get_windowed_coverage_and_diagnostic_reads.sh
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1677557222"
    memory: "4 GiB"
    preemptible: 3
  }
  output {
    String message = read_string(stdout())
  }
}


