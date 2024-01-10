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
    String sample_id
    String output_dir
    Int interval
    Int window_size
    Int min_qual
    #coverage summary parameters
    Float accessibility_threshold
    Float mapq_threshold
    File accessibility_mask_file
    File mapq_file
    File gc_content_file
    String sample_group_id
    String species
    String runtime_zones
  }

  call WindowedCoverage {
    input:
      input_bam = input_bam,
      input_bai = input_bai,
      sample_id = sample_id,
      output_dir = output_dir,
      interval = interval,
      window_size = window_size,
      min_qual = min_qual,
      runtime_zones = runtime_zones
  }

  call CoverageSummary {
    input:
      coverage_tarball = WindowedCoverage.output_gz,
      accessibility_threshold = accessibility_threshold,
      mapq_threshold = mapq_threshold,
      accessibility_mask_file = accessibility_mask_file,
      mapq_file = mapq_file,
      gc_content_file = gc_content_file,
      sample_group_id = sample_group_id,
      output_dir = output_dir,
      sample_id = sample_id,
      runtime_zones = runtime_zones
  }

  call CoverageHMM {
    input:
      coverage_tarball = CoverageSummary.output_gz,
      gc_content_file = gc_content_file,
      coverage_gc = CoverageSummary.coverage_output,
      coverage_variance = CoverageSummary.variance_output,
      mapq_file = mapq_file,
      mapq_threshold = mapq_threshold,
      species =  species,
      output_dir = output_dir,
      sample_id = sample_id,
      runtime_zones = runtime_zones
  }
  
  output {
    File output_gz = CoverageHMM.output_gz
    File coverage_variance = CoverageSummary.variance_output
  }
}

task WindowedCoverage {
  input {
    File input_bam
    File input_bai
    String sample_id
    String output_dir
    Int interval
    Int window_size
    Int min_qual

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    Int disk = ceil(size(input_bam, "GiB")) + 50
    Int preemptible = 3
    String runtime_zones
  }
  meta {
    description: "Compute aligned read counts across the genome in 300 bp windows. The output is a set of Windowed count reads file, 1 row per window, 1 file per sample - compressed in a tar.gz file to keep the directory structure."
  }
  parameter_meta {
    input_bam: "The input BAM file"
    input_bai: "The BAM file's corresponding index file"
    sample_id: "The sample name. This is typically the same as the input bam file name but without the extension."
    output_dir: "The output directory name. This is set to `coverage` by default."
    interval: "The interval size (bp) to use for the coverage calculations. Default is 300."
    window_size: "The window size (bp) to use for the coverage calculations. Default is 300."
    min_qual: "The minimum quality to use for the coverage calculations. Default is 20."
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    runtime_zones: "The ordered list of zone preference for running in GCP"
  }
  command <<<
    set -x
    echo "Calculating coverage for: " 
    basename ~{input_bam}
    echo "Sample Name: ~{sample_id}"
    echo "Current directory: " 
    pwd
    ls -lht
    # Make the output directory
    mkdir -p ~{output_dir}
    # Start the conda env
    source activate cnv37 
    # Get the coverage data for each chromosome
    allchrom=(2L 2R 3L 3R X)
    for chrom in ${allchrom[@]}
    do
      #create the directory structure needed (this was likely pre-provisioned in the baremetal version of this pipeline)
      mkdir -p ~{output_dir}/${chrom}/coveragelogs/
      python /cnv/scripts/counts_for_HMM.py \
            ~{input_bam} \
            $chrom \
            ~{interval} ~{window_size} ~{min_qual} \
            ~{output_dir}/${chrom}/counts_for_HMM_~{sample_id}_${chrom}_output.csv \
            > ~{output_dir}/${chrom}/coveragelogs/counts_for_HMM_~{sample_id}_${chrom}.log 2>&1
    done
    # Compress the output directory
    tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
    ls -lht
  >>>
  runtime {
    docker: docker
    memory: ram
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
    zones: runtime_zones
  }
  output {
    #Compressed output directory
    File output_gz = "~{output_dir}.tar.gz"
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
    gc_content_file: "The gc content file to use for the coverage summary calculations"
    sample_group_id: "The sample group id to use for the coverage summary calculations"
    sample_id: "The sample name. This is typically the same as the input bam file name but without the extension."
    output_dir: "The output directory name. This is set to `coverage` by default."
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    runtime_zones: "The ordered list of zone preference for running in GCP"
  }
  input {
    File coverage_tarball
    Float accessibility_threshold
    Float mapq_threshold
    File accessibility_mask_file
    File mapq_file
    File gc_content_file
    String sample_group_id
    String output_dir
    String sample_id
    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    # TODO: Make disk space dynamic based on input size
    Int disk = 70
    Int preemptible = 3
    String runtime_zones
  }
  String acc_threshold = sub(accessibility_threshold + "_", "[.]", "")
  String m_threshold = sub(mapq_threshold + "_", "[.]", "")
  String coverage_output_filename = "median_coverage_by_GC_masked_" + acc_threshold + m_threshold + sample_group_id + ".csv"
  String variance_output_filename = "coverage_variance_masked_" + acc_threshold + m_threshold + sample_group_id + ".csv"
  command <<<
    set -x
    echo "Current directory: " 
    pwd
    #unzip the tarball
    tar -zxvf ~{coverage_tarball}
    #cd /cnv/scripts/
    echo "Creating temporary manifest: "
    echo ~{sample_id} > manifest.txt
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
            manifest.txt \
            ~{gc_content_file} \
            ~{output_dir} \
            ~{sample_group_id} \
            > calculate_mean_coverage_by_GC_~{sample_group_id}.log 2>&1
    ls -lht
    tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
    echo "Numbers: ~{acc_threshold} ~{m_threshold} ~{sample_group_id}"
  
  >>>
  runtime {
    docker: docker
    memory: ram
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
    zones: runtime_zones
  }
  output {
    File logs = "calculate_mean_coverage_by_GC_~{sample_group_id}.log"
    File coverage_output = "~{output_dir}/~{coverage_output_filename}"
    File variance_output = "~{output_dir}/~{variance_output_filename}"
    File output_gz = "~{output_dir}.tar.gz"
  }
}

task CoverageHMM {
  meta {
    description: "Runs an mHMM on coverage counts obtained from a bamfile, using normalisation based on the mean coverage per GC bin over the whole genome of the individual. The outputs are normalized coverage values, 1 per window, per sample. Copy number state, 1 per window, per sample."
  }
  parameter_meta {
    coverage_tarball: "The input tarball containing the coverage files"
    gc_content_file: "The gc content file to use for the HMM calculations"
    coverage_gc: "The coverage gc file to use for the HMM calculations"
    coverage_variance: "The coverage variance file to use for the HMM calculations"
    mapq_file: "The mapq file to use for the HMM calculations"
    mapq_threshold: "The mapq threshold to use for the HMM calculations. Default is 0.5."
    species: "The species of the samples specified in the manifest."
    sample_id: "The sample name. This is typically the same as the input bam file name but without the extension."
    docker: "(optional) the docker image containing the runtime environment for this task"
    ram: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    runtime_zones: "The ordered list of zone preference for running in GCP"
  }
  input {
    File coverage_tarball
    File gc_content_file
    File coverage_gc
    File coverage_variance
    File mapq_file
    Float mapq_threshold
    String species
    String output_dir
    String sample_id
    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
    String ram = "8000 MiB"
    Int cpu = 16
    Int disk = 70
    Int preemptible = 3
    String runtime_zones
  }
  command <<<
    set -x
    echo "Current directory: " 
    pwd
    #unzip the tarball
    tar -zxvf ~{coverage_tarball}
    echo "Creating temporary manifest: "
    echo ~{sample_id} > manifest.txt
    ls -lht
    allchrom=(2L 2R 3L 3R X)
    # Activate the conda environment
    source activate cnv37 
    for chrom in ${allchrom[@]}
    do
      mkdir -p ~{output_dir}/$chrom/HMM_logs_~{species}/
      mkdir -p ~{output_dir}/$chrom/HMM_output/
      # run an mHMM on coverage counts
      python /cnv/scripts/HMM_process.py \
        manifest.txt \
        $chrom \
        ~{output_dir} \
        ~{gc_content_file} \
        ~{coverage_gc} \
        ~{coverage_variance} \
        ~{mapq_file} \
        ~{mapq_threshold} \
        > ~{output_dir}/$chrom/HMM_logs_~{species}/HMM_$chrom.log 2>&1
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
    zones: runtime_zones
  }
  output {
    #log outputs
    Array[File] logs = glob("*.log")
    File output_gz = "~{output_dir}.tar.gz"
  }
}
