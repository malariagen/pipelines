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
    File output_gz = WindowedCoverage.output_gz
  }
}

task WindowedCoverage {
  input {
    File input_bam
    String sample_name
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
    echo "Calculating coverage for: " 
    basename ~{input_bam}
    echo "Current directory: " 
    pwd
    ls -lht
    cd /cnv/scripts/
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
      #mkdir -p ${coveragefolder}/${chrom}/coveragelogs/
      python /cnv/scripts/counts_for_HMM.py \
            ~{input_bam} \
            $chrom \
            300 300 10 \
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


