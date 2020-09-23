version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task SelectVariants {

  input {
    File sample_vcf
    File alleles_vcf
    String output_basename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }
#TODO: is it better to pass individual files or combine into one alleles vcf?
#TODO: do we need to index here or will the file be available? (needs an index file for the sample vcf, but not the intervals vcf)
  command {
    # generate index file if needed
    gatk IndexFeatureFile \
      --input ~{sample_vcf}

    # subset the VCF
    gatk SelectVariants \
      -R REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      -V ~{sample_vcf} \
      --select-type-to-include SNP \
      -O ~{output_basename}.vcf \
      --remove-unused-alternates true \
      --intervals ~{alleles_vcf} \
      --intervals ~{alleles_vcf}
  }

  runtime {
    docker: runTimeSettings.gatk_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
    File sample_subset_vcf = "~{output_basename}.subset.vcf"
  }
}

# TODO: try running this
# TODO: add  --max-coverage=@@TODO
task WhatsHapPhase {
  input {
    File sample_bam
    File sample_subset_vcf
    String output_basename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    whatshap phase \
      -o ~{output_basename}_phased.vcf \
      --reference ~{reference.ref_fasta} \
      ~{sample_subset_vcf} ~{sample_bam}
  }

  runtime {
    docker: runTimeSettings.whatshap_docker # us.gcr.io/broad-gotc-prod/malariagen/whatshap:0.0
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
    File sample_phased_vcf = "~{output_basename}.phased.vcf"
  }
}

# TODO: what is chromosome lengths??
# TODO: run this. What is the output??
task WhatsHapStats {
  input {
   File sample_phased_vcf
   String output_basename
   RunTimeSettings runTimeSettings
  }

  command {
    whatshap stats \
      --chr-lengths=CHR_LENGTHS \
      --tsv=~{output_basename}.stats.tsv \
      --gtf=~{output_basename}.blocks.gtf \
      ~{sample_phased_vcf}
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