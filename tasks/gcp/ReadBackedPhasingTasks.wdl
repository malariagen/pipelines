version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task SelectVariants {

  input {
    File input_vcf
    Array[File] phased_sites_vcfs
    Array[File] phased_sites_indicies
    String output_basename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(input_vcf, "GiB") * 3)

  command {
    # generate index file if needed
    gatk IndexFeatureFile \
      --input ~{input_vcf}

    # subset the VCF
    gatk SelectVariants \
      -R ~{reference.ref_fasta} \
      -V ~{input_vcf} \
      --select-type-to-include SNP \
      -O ~{output_basename}.subset.vcf \
      --remove-unused-alternates true \
      --intervals ~{sep=' --intervals ' phased_sites_vcfs}
  }

  runtime {
    docker: runTimeSettings.gatk_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "7.5 GiB"
    disks: "local-disk ${disk_size} HDD"
  }

  output {
    File subset_vcf = "~{output_basename}.subset.vcf"
  }
}

# TODO: try running this
# TODO: add  --max-coverage=@@TODO
# TODO: update docker tag to be tool version number
task WhatsHapPhase {
  input {
    File input_bam
    File input_bam_index
    File subset_vcf
    String output_basename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    whatshap phase \
      -o ~{output_basename}_phased.vcf \
      --reference ~{reference.ref_fasta} \
      ~{subset_vcf} ~{input_bam}
  }

  runtime {
    docker: runTimeSettings.whatshap_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
    File phased_vcf = "~{output_basename}.phased.vcf"
  }
}

# TODO: what is chromosome lengths??
# TODO: run this. What is the output??
task WhatsHapStats {
  input {
   File phased_vcf
   String output_basename
   RunTimeSettings runTimeSettings
  }

  command {
    whatshap stats \
      --chr-lengths=CHR_LENGTHS \
      --tsv=~{output_basename}.stats.tsv \
      --gtf=~{output_basename}.blocks.gtf \
      ~{phased_vcf}
  }

  runtime {
    docker: runTimeSettings.whatshap_docker # us.gcr.io/broad-gotc-prod/malariagen/whatshap:0.0
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
    File whats_hap_stats = "~{output_basename}.whatshap_stats"
  }
}