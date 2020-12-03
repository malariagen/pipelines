version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"

task BgzipAndTabix {
    input {
        File input_vcf
        String output_basename
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
        docker: runTimeSettings.bcftools_docker
        preemptible: runTimeSettings.preemptible_tries
        cpu: "1"
        memory: "3.75 GiB"
    }
    output {
        File vcf = "~{output_basename}.vcf.gz"
        File vcf_index = "~{output_basename}.vcf.gz.tbi"
    }
}