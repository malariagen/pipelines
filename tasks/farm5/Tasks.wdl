version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"

task BgzipAndTabix {
  input {
    File input_vcf
    String output_basename

    String singularity_image = "bcftools.1.11.sif"
    Int num_cpu = 1
    Int memory = 32000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    # note that bgzip has an option (-i) to index the bgzipped output, but this file is not a tabix file
    # note also that we use '-c' so that bgzip doesn't create the bgzipped file in place, rather it's in a location
    # where it's easy to output from the task.
    bgzip -c ~{input_vcf} > ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File vcf = "~{output_basename}.vcf.gz"
    File vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task Tabix {
  input {
    File input_file

    String singularity_image = "bcftools.1.11.sif"
    Int num_cpu = 1
    Int memory = 32000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  String local_file = basename(input_file)

  command {
    # Localize the passed input_file to the working directory so when the
    # newly created index file doesn't get delocalized with the long path.
    cp ~{input_file} ~{local_file}
    tabix ~{local_file}
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    # output the path to the copied local file AND the created index so they are side by side.
    File output_file = local_file
    File output_index_file = "~{local_file}.tbi"
  }
}