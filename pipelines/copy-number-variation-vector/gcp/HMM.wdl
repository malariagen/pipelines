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
    Int interval
    Int window_size
    Int min_qual
  }
  call WindowedCoverage {
    input:
      input_bam = input_bam,
      sample_name = sample_name,
      output_dir = output_dir,
      interval = interval,
      window_size = window_size,
      min_qual = min_qual
  }

  call CoverageSummary {
    input:
      coverage_tarball = WindowedCoverage.output_gz,
      accessibility_threshold = 0.9
  }
  
  output {
    File output_gz = WindowedCoverage.output_gz
  }
}

task WindowedCoverage {
  input {
    File input_bam
    String sample_name
    String output_dir
    Int interval
    Int window_size
    Int min_qual

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    # TODO: Make disk space dynamic based on input size
    Int disk = 70
    Int preemptible = 3
  }
  meta {
    description: "Compute aligned read counts across the genome in 300 bp windows. The output is a set of Windowed count reads file, 1 row per window, 1 file per sample - compressed in a tar.gz file to keep the directory structure."
  }
  parameter_meta {
    input_bam: "The input bam file"
    sample_name: "The sample name. This is typically the same as the input bam file name but without the extension."
    output_dir: "The output directory name. This is set to `coverage` by default."
    interval: "The interval size (bp) to use for the coverage calculations. Default is 300."
    window_size: "The window size (bp) to use for the coverage calculations. Default is 300."
    min_qual: "The minimum quality to use for the coverage calculations. Default is 20."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  command <<<
    set -x
    echo "Calculating coverage for: " 
    basename ~{input_bam}
    echo "Current directory: " 
    pwd
    ls -lht
    #cd /cnv/scripts/
    ls -lht /cnv
    #bash get_windowed_coverage_and_diagnostic_reads.sh ~{input_bam} ~{sample_name} ~{output_dir}
    mkdir -p ~{output_dir}
    # Start the conda env
    source activate cnv37 
    # Get the coverage data
    allchrom=(2L 2R 3L 3R X)
    #coveragefolder=~{output_dir}/coverage
    for chrom in ${allchrom[@]}
    do
      #create the directory structure needed (this was likely pre-provisioned in the baremetal version of this pipeline)
      mkdir -p ~{output_dir}/${chrom}/coveragelogs/
      python /cnv/scripts/counts_for_HMM.py \
            ~{input_bam} \
            $chrom \
            ~{interval} ~{window_size} ~{min_qual} \
            ~{output_dir}/${chrom}/counts_for_HMM_${samplename}_${chrom}_output.csv \
            > ~{output_dir}/${chrom}/coveragelogs/counts_for_HMM_${samplename}_${chrom}.log 2>&1
    done
    tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
    ls -lht
  >>>
  runtime {
    docker: docker
    memory: ram
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    #Compressed output directory
    File output_gz = "~{output_dir}.tar.gz"
    #String message = read_string(stdout())
    #Return an array of Files representing the contents of the output directory
    # Array[File] coverage_files = glob("~{output_dir}/*.csv")
    #Return an array of Files representing the logs of the coverage calculations
    # Array[File] coverage_logs = glob("~{output_dir}/*.log")
  }
}

task CoverageSummary {
  meta {
    description: "Analyzes the read count for each sample. Computes a GC normalization table that describes the mean read count for each percentile of GC. Computes overall coverage variance. The output is a summary statistics file, 1 per sample"
  }
  input {
    File coverage_tarball
    Int accessibility_threshold

    # Accessibility Threshold = 0.9 (default, hard-coded in the call to calculate_median_coverage_by_GC)
    # Mapq Threshold = 0.5 (default, hard-coded in the call to calculate_median_coverage_by_GC)
    # Accessibility Mask File
    # Mapq File
    # Sample Manifest
    # GC Content File (GC composition for the Agam genom)
    # Sample Group ID (Appears to be a string only used to append to the output filename)
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
    echo "Calculating stats for: " 
    #basename ~{input_bam}
    echo "Current directory: " 
    pwd
    #unzip the tarball
    tar -zxvf ~{coverage_tarball}
    #cd /cnv/scripts/
    ls -lht
    
  >>>
}
