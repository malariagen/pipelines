version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task UnifiedGenotyper {
  input {
    File input_bam
    File input_bam_index
    File alleles_vcf
    File alleles_vcf_index
    String output_vcf_filename

    String docker_tag = "malaria-gatk3@sha256:11dcafb2c5b574c8313942874ec8449b536adcc37c00bad149f1ef1a45012a28"
    Int num_cpu = 4
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx~{memory}m \
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
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String docker_tag = "malaria-samplevcftozarr@sha256:1baec1f2b253311bf834f8c5bf8c8169e765f738ac1972100027cbfa28329f2f"
    Int num_cpu = 2
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

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
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_log_file = output_log_file_name
    File zarr_output = "~{output_zarr_file_name}.zip"
  }
}

