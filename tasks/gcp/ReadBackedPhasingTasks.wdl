version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task SelectVariants {

  input {
    File sample_zarr
    File called_sites_zarr
    File phased_sites_zarr
    String output_basename
    String contig
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int mem_size = ceil(size(sample_zarr, "GiB") * 3)
  Int disk_size = ceil(size(sample_zarr, "GiB") * 2)

  command {
    python /tools/sample_select_variants.py \
      --sample-genotypes ~{sample_zarr} \
      --sites-called ~{called_sites_zarr} \
      --sites-selected ~{phased_sites_zarr} \
      --output ~{output_basename}.subset.vcf \
      --contig ~{contig} \
      --progress
  }

  runtime {
    docker: runTimeSettings.select_variants_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
    disks: "local-disk 100 HDD"
  }

  output {
    File subset_vcf = "~{output_basename}.subset.vcf"
  }
}
    #File subset_vcf_index = "~{output_basename}.subset.vcf.gz.tbi"


# TODO: add  --max-coverage=@@TODO
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
      -o ~{output_basename}.phased.vcf \
      --reference ~{reference.ref_fasta} \
      ~{subset_vcf} ~{input_bam}
  }

  runtime {
    docker: runTimeSettings.whatshap_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "4"
    memory: "30 GiB"
  }

  output {
    File phased_vcf = "~{output_basename}.phased.vcf"
  }
}

# TODO: what is chromosome lengths??
task WhatsHapStats {
  input {
   File phased_vcf
   String output_basename
   RunTimeSettings runTimeSettings
  }
      #--chr-lengths=CHR_LENGTHS \
  command {
    whatshap stats \
      ~{phased_vcf} \
      --tsv=~{output_basename}.stats.tsv \
      --gtf=~{output_basename}.blocks.gtf
  }

  runtime {
    docker: runTimeSettings.whatshap_docker # us.gcr.io/broad-gotc-prod/malariagen/whatshap:0.0
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
    File whats_hap_stats_tsv = "~{output_basename}.stats.tsv"
    File whats_hap_blocks_gtf = "~{output_basename}.blocks.gtf"
  }
}