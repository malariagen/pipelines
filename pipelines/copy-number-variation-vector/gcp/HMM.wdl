version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2021
##
## This WDL pipeline implements the Read-Backed component of the Mospquito Phasing Pipeline as described in
## https://github.com/malariagen/pipelines/blob/master/docs/specs/phasing-vector.md
##


workflow HMM {
  String pipeline_version = "1.0.0"

  input {
    #windowed coverage parameters
    File input_bam
    File input_bai
    String sample_name
    String output_dir
    Int interval
    Int window_size
    Int min_qual
    #coverage summary parameters
    Float accessibility_threshold
    Float mapq_threshold
    File accessibility_mask_file
    File mapq_file
    File sample_manifest
    File gc_content_file
    String sample_group_id
  }

  call WindowedCoverage {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      sample_name = sample_name,
      output_dir = output_dir,
      interval = interval,
      window_size = window_size,
      min_qual = min_qual
  }

  call CoverageSummary {
    input:
      coverage_tarball = WindowedCoverage.output_gz,
      accessibility_threshold = accessibility_threshold,
      mapq_threshold = mapq_threshold,
      accessibility_mask_file = accessibility_mask_file,
      mapq_file = mapq_file,
      sample_manifest = sample_manifest,
      gc_content_file = gc_content_file,
      sample_group_id = sample_group_id,
      output_dir = output_dir
  }
  
  output {
    File output_gz = CoverageSummary.output_gz
  }
}

task WindowedCoverage {
  input {
    File input_bam
    File input_bai
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
    input_bam: "The input BAM file"
    input_bai: "The BAM file's corresponding index file"
    sample_name: "The sample name. This is typically the same as the input bam file name but without the extension."
    output_dir: "The output directory name. This is set to `coverage` by default."
    interval: "The interval size (bp) to use for the coverage calculations. Default is 300."
    window_size: "The window size (bp) to use for the coverage calculations. Default is 300."
    min_qual: "The minimum quality to use for the coverage calculations. Default is 20."
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  command <<<
    set -x
    echo "Calculating coverage for: " 
    basename ~{input_bam}
    echo "Sample Name: ~{sample_name}"
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
            ~{output_dir}/${chrom}/counts_for_HMM_~{sample_name}_${chrom}_output.csv \
            > ~{output_dir}/${chrom}/coveragelogs/counts_for_HMM_~{sample_name}_${chrom}.log 2>&1
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
  parameter_meta {
    coverage_tarball: "The input tarball containing the coverage files"
    accessibility_threshold: "The accessibility threshold to use for the coverage summary calculations. Default is 0.9."
    mapq_threshold: "The mapq threshold to use for the coverage summary calculations. Default is 0.5."
    accessibility_mask_file: "The accessibility mask file to use for the coverage summary calculations"
    mapq_file: "The mapq file to use for the coverage summary calculations"
    sample_manifest: "The sample manifest file to use for the coverage summary calculations. This is simply a list of all the sample names."
    gc_content_file: "The gc content file to use for the coverage summary calculations"
    sample_group_id: "The sample group id to use for the coverage summary calculations"
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  input {
    File coverage_tarball
    Float accessibility_threshold
    Float mapq_threshold
    File accessibility_mask_file
    File mapq_file
    File sample_manifest
    File gc_content_file
    String sample_group_id
    String output_dir
    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    # TODO: Make disk space dynamic based on input size
    Int disk = 70
    Int preemptible = 3
    String coverage_output_filename = ""
  }
  command <<<
    set -x
    echo "Current directory: " 
    pwd
    #unzip the tarball
    tar -zxvf ~{coverage_tarball}
    #cd /cnv/scripts/
    ls -lht
    allchrom=(2L 2R 3L 3R X)
    source activate cnv37 
    # Calculate median coverage by GC 
    echo "Calculating median and variance of coverage by GC bin for sample group ~{sample_group_id}"
    python /cnv/scripts/calculate_median_coverage_by_GC.py \
            ~{accessibility_threshold} \
            ~{accessibility_mask_file} \
            ~{mapq_threshold} \
            ~{mapq_file} \
            ~{sample_manifest} \
            ~{gc_content_file} \
            ~{output_dir} \
            ~{sample_group_id} \
            > calculate_mean_coverage_by_GC_~{sample_group_id}.log 2>&1
    ls -lht
    tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
    #output_filename = working_folder + '/median_coverage_by_GC_masked_' + sub('\.', '', str(accessibility_threshold)) + '_' + sub('\.', '', str(mapq0_threshold)) + '_' + output_file_key + '.csv'
    #output_variance_filename = working_folder + '/coverage_variance_masked_' + sub('\.', '', str(accessibility_threshold)) + '_' + sub('\.', '', str(mapq0_threshold)) + '_' + output_file_key + '.csv'
    #TODO: Create output variables for the output of this script
    acc_threshold=~{accessibility_threshold}
    m_threshold=~{mapq_threshold}
    sg_id=~{sample_group_id}
    coverage_output_str="~{output_dir}/median_coverage_by_GC_masked_${acc_threshold//./_}_${m_threshold//./_}_${sg_id//./_}.csv"
    echo "Coverage output string resolved: ${coverage_output_str}"
    coverage_output_filename="${coverage_output_str}"
  >>>
  runtime {
    docker: docker
    memory: ram
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    #log output
    #${orig//[xyz]/_} ${~{accessibility_threshold}//./} ${~{mapq_threshold}//./}
    #select_first(glob("~{output_dir}/median_coverage_by_GC_masked_*_~{sample_group_id}.csv"))
    File logs = "calculate_mean_coverage_by_GC_~{sample_group_id}.log"
    File coverage_output = "~{coverage_output_filename}"
    #File variance_output = "~{output_dir}/coverage_variance_masked_*_~{sample_group_id}.csv"
    File output_gz = "~{output_dir}.tar.gz"
  }
}

task CoverageHMM {
  meta {
    description: "Runs an mHMM on coverage counts obtained from a bamfile, using normalisation based on the mean coverage per GC bin over the whole genome of the individual. The outputs are normalized coverage values, 1 per window, per sample. Copy number state, 1 per window, per sample."
  }
  parameter_meta {
    sample_manifest: "The sample manifest file to use for the HMM calculations. This is simply a list of all the sample names."
    chrom: "The chromosome to run the coverage HMM on"
    coverage_tarball: "The input tarball containing the coverage files"
    gc_content_file: "The gc content file to use for the HMM calculations"
    coverage_gc: "The coverage gc file to use for the HMM calculations"
    coverage_variance: "The coverage variance file to use for the HMM calculations"
    mapq_file: "The mapq file to use for the HMM calculations"
    mapq_threshold: "The mapq threshold to use for the HMM calculations. Default is 0.5."
    species: "The species of the samples specified in the manifest."
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }
  input {
    File sample_manifest
    String chrom
    File coverage_tarball
    File gc_content_file
    File coverage_gc
    File coverage_variance
    File mapq_file
    Float mapq_threshold
    String species

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
    #unzip the tarball
    tar -zxvf ~{coverage_tarball}
    #cd /cnv/scripts/
    ls -lht
    allchrom=(2L 2R 3L 3R X)
    # Activate the conda environment
    source activate cnv37 
    # Calculate median coverage by GC 
    python /cnv/scripts/HMM_process.py \
      ~{sample_manifest} \
      ~{chrom} \
      coverage \
      ~{gc_content_file} \
      ~{coverage_gc} \
      ~{coverage_variance} \
      ~{mapq_file} \
      ~{mapq_threshold} \
      > coverage/~{chrom}/HMM_logs_~{species}/HMM_~{chrom}.log 2>&1

      #  $manifest \
      #  $chrom \
      #  $coveragefolder \
      #  $GC_content_file \
      #  $coverage_by_GC_file \
      #  $coverage_variance_file \
      #  $mapq_prop_file \
      #  0.5 \
      #  > ${coveragefolder}/${chrom}/HMM_logs_${species}/HMM_${chrom}.log 2>&1
  >>>
}
