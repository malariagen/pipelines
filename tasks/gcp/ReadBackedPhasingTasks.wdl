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

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/sampleselectvariants:1.0"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
    String runtime_zones
  }

  Int mem_size = ceil(size(sample_zarr, "GiB") * 3)
  Int disk_size = ceil(size(sample_zarr, "GiB") + size(called_sites_zarr, "GiB")+ size(phased_sites_zarr, "GiB")) * 2 + 20

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
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }

  output {
    File subset_vcf = "~{output_basename}.subset.vcf"
  }
}


task WhatsHapPhase {
  input {
    File input_bam
    File input_bam_index
    File subset_vcf
    File subset_vcf_index
    String output_filename
    String? contig
    Int internal_downsampling = 15
    ReferenceSequence reference

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/whatshap:1.0"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 2
    RunTimeSettings runTimeSettings
    String runtime_zones
  }

  Int disk_size = ceil(size(subset_vcf, "GiB") + size(subset_vcf_index, "GiB") + size(input_bam, "GiB") + size(input_bam_index, "GiB")) * 5 + 20

  command {
    touch ~{input_bam_index}
    touch ~{subset_vcf_index}
    whatshap phase \
      -o ~{output_filename} \
      --reference ~{reference.ref_fasta} \
      ~{"--chromosome " + contig} \
      --internal-downsampling=~{internal_downsampling} \
      ~{subset_vcf} ~{input_bam}
  }

  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "50 GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }

  output {
    File phased_vcf = output_filename
  }
}

task WhatsHapStats {
  input {
    File phased_vcf
    File phased_vcf_index
    String output_basename
    ReferenceSequence reference

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/whatshap:1.0"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 2
    RunTimeSettings runTimeSettings
    String runtime_zones
  }
  # TODO - figure out how to make proper use of OPTIONAL reference.ref_chr_lengths
  command {
    touch ~{phased_vcf_index}
    whatshap stats \
      ~{phased_vcf} \
      --chr-lengths=~{reference.ref_chr_lengths} \
      --tsv=~{output_basename}.stats.tsv \
      --gtf=~{output_basename}.blocks.gtf
  }

  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "7.5 GiB"
    zones: runtime_zones
  }

  output {
    File whats_hap_stats_tsv = "~{output_basename}.stats.tsv"
    File whats_hap_blocks_gtf = "~{output_basename}.blocks.gtf"
  }
}
