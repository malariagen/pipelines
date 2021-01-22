version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task UnifiedGenotyper {
  input {
    File input_bam
    File input_bam_index
    File alleles_vcf
    File alleles_vcf_index
    String output_vcf_filename

    String docker_tag = "broadinstitute/gatk3:3.7-0"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 4
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 4) + 20

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T UnifiedGenotyper \
          -I ~{input_bam} \
          --alleles ~{alleles_vcf} \
          -R ~{reference.ref_fasta} \
          --out ~{output_vcf_filename} \
          --genotype_likelihoods_model BOTH \
          --genotyping_mode GENOTYPE_GIVEN_ALLELES \
          --heterozygosity 0.015 \
          --heterozygosity_stdev 0.05 \
          --indel_heterozygosity 0.001 \
          --downsampling_type BY_SAMPLE \
          -dcov 250 \
          --output_mode EMIT_ALL_SITES \
          --min_base_quality_score 17 \
          -stand_call_conf 0.0 \
          -contamination 0.0 \
          -A DepthPerAlleleBySample \
          -A RMSMappingQuality \
          -XA Coverage \
          -XA ExcessHet \
          -XA InbreedingCoeff \
          -XA MappingQualityZero \
          -XA HaplotypeScore \
          -XA SpanningDeletions \
          -XA FisherStrand \
          -XA StrandOddsRatio \
          -XA ChromosomeCounts \
          -XA BaseQualityRankSumTest \
          -XA MappingQualityRankSumTest \
          -XA QualByDepth \
          -XA ReadPosRankSumTest
  }
  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "15 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = "~{output_vcf_filename}.idx"
  }
}

task VcfToZarr {
  input {
    File input_vcf
    String sample_id
    String output_zarr_file_name
    String output_log_file_name

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/samplevcftozarr:1.1"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_vcf, "GiB")) * 4) + 20

  command {
    python /tools/sample_vcf_to_zarr.py \
        --input ~{input_vcf} \
        --output ~{output_zarr_file_name} \
        --sample ~{sample_id} \
        --field variants/MQ \
        --field calldata/GT \
        --field calldata/GQ \
        --field calldata/AD \
        --contig 2R \
        --contig 2L \
        --contig 3R \
        --contig 3L \
        --contig X \
        --contig Y_unplaced \
        --contig UNKN \
        --log ~{output_log_file_name} \
        --zip
  }
  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "7.5 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_log_file = output_log_file_name
    File zarr_output = "~{output_zarr_file_name}.zip"
  }
}

