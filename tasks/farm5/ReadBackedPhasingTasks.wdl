version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task SelectVariants {

  input {
    File sample_zarr
    File called_sites_zarr
    File phased_sites_zarr
    String output_basename
    String contig

    String singularity_image = "sampleselectvariants.1.0.sif"
    Int num_cpu = 1
    Int memory = 5000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = "whatshap.1.0.sif"
    Int num_cpu = 2
    Int memory = 32000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = "whatshap.1.0.sif"
    Int num_cpu = 2
    Int memory = 9000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
    File whats_hap_stats_tsv = "~{output_basename}.stats.tsv"
    File whats_hap_blocks_gtf = "~{output_basename}.blocks.gtf"
  }
}
