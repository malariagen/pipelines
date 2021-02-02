version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"

task BgzipAndTabix {
  input {
    File input_vcf
    String output_basename

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/bcftools:1.11"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(input_vcf, "GiB")) * 3 + 20

  command {
    # note that bgzip has an option (-i) to index the bgzipped output, but this file is not a tabix file
    # note also that we use '-c' so that bgzip doesn't create the bgzipped file in place, rather it's in a location
    # where it's easy to output from the task.
    bgzip -c ~{input_vcf} > ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  }
  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File vcf = "~{output_basename}.vcf.gz"
    File vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task Tabix {
  input {
    File input_file

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/bcftools:1.11"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(input_file, "GiB")) * 2 + 20
  String local_file = basename(input_file)

  command {
    # Localize the passed input_file to the working directory so when the
    # newly created index file doesn't get delocalized with the long path.
    cp ~{input_file} ~{local_file}
    tabix ~{local_file}
  }
  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    # output the path to the copied local file AND the created index so they are side by side.
    File output_file = local_file
    File output_index_file = "~{local_file}.tbi"
  }
}