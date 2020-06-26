version 1.0

import "../../structs/ReferenceSequence.wdl"

workflow AmpliconSNPCallingParasite {
  input {
    File input_bam
    File input_bam_index
    String output_basename
    ReferenceSequence reference

    File mpileup_targets_vcf
    Int? mpileup_min_base_quality
    Int? mpileup_max_depth

    Int? call_ploidy

    Int? filter_min_depth
    Int? filter_min_quality
    Int? filter_min_mapping_quality

    Int? preemptible_tries
    String docker
  }

  call CallAndFilterSNPsFromPileups {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_basename = output_basename,
      reference = reference,
      mpileup_targets_vcf = mpileup_targets_vcf,
      mpileup_min_base_quality = mpileup_min_base_quality,
      mpileup_max_depth = mpileup_max_depth,
      call_ploidy = call_ploidy,
      filter_min_depth = filter_min_depth,
      filter_min_quality = filter_min_quality,
      filter_min_mapping_quality = filter_min_mapping_quality,
      preemptible_tries = preemptible_tries,
      docker = docker
  }

  output {
    File output_pileup_vcf = CallAndFilterSNPsFromPileups.output_pileup_vcf
    File output_unfiltered_calls_vcf = CallAndFilterSNPsFromPileups.output_unfiltered_calls_vcf
    File output_filtered_calls_vcf = CallAndFilterSNPsFromPileups.output_filtered_calls_vcf
  }
}

task CallAndFilterSNPsFromPileups {
  input {
    File input_bam
    File input_bam_index
    String output_basename
    ReferenceSequence reference

    File mpileup_targets_vcf
    Int? mpileup_min_base_quality
    Int? mpileup_max_depth

    Int? call_ploidy

    Int? filter_min_depth
    Int? filter_min_quality
    Int? filter_min_mapping_quality

    Int? preemptible_tries
    String docker
  }

  Float bam_size = size(input_bam, "GiB")
  Float ref_size = size(reference.ref_fasta, "GiB") + size(reference.ref_fasta_index, "GiB") + size(reference.ref_dict, "GiB")
  Float disk_multiplier = 2.5
  Int disk_size = ceil(bam_size + ref_size + (disk_multiplier * bam_size) + 20)

  String output_pileup_vcf_ = output_basename + ".pileup.vcf"
  String output_unfiltered_calls_vcf_ = output_basename + ".unfiltered.calls.vcf"
  String output_filtered_calls_vcf_ = output_basename + ".filtered.calls.vcf"

  command {
    set -euo pipefail

    echo "* * * * ~{default="2" call_ploidy}" > ploidy.txt

    bcftools mpileup \
      --min-BQ ~{default="20" mpileup_min_base_quality} \
      --annotate FORMAT/AD,FORMAT/DP \
      --max-depth ~{default="50000" mpileup_max_depth} \
      --targets-file ~{mpileup_targets_vcf} \
      --fasta-ref ~{reference.ref_fasta} \
      --output-type v \
      ~{input_bam} \
      > ~{output_pileup_vcf_}

    bcftools call \
      --no-version \
      --multiallelic-caller \
      --keep-alts \
      --skip-variants indels \
      --ploidy-file ploidy.txt \
      --output-type u \
      < ~{output_pileup_vcf_} \
      > ~{output_unfiltered_calls_vcf_}

    bcftools filter \
      --mode + \
      --soft-filter LowDepth \
      --exclude FORMAT/DP\<~{default="8" filter_min_depth} \
      --output-type v \
      < ~{output_unfiltered_calls_vcf_} |
    bcftools filter \
      --mode + \
      --soft-filter LowQual \
      --exclude '%QUAL<~{default="15" filter_min_quality} || MQ<~{default="20" filter_min_mapping_quality}' \
      --output-type v \
      > ~{output_filtered_calls_vcf_}
  }

  runtime {
    docker: docker
    preemptible: select_first([preemptible_tries, 5])
    cpu: "1"
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_pileup_vcf = "~{output_pileup_vcf_}"
    File output_unfiltered_calls_vcf = "~{output_unfiltered_calls_vcf_}"
    File output_filtered_calls_vcf = "~{output_filtered_calls_vcf_}"
  }
}
