version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {
    Array[File] phased_sample_vcfs
    Array[File] phased_sample_vcf_indices
    String project_id

    String singularity_image = "bcftools.1.11.sif"
    Int num_cpu = 1
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  #TODO - touch all the index files to avoid the annoying warning
  command {
    bcftools merge \
      -o ~{project_id}_merged.vcf \
      ~{sep=' ' phased_sample_vcfs}
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File merged_vcf = "~{project_id}_merged.vcf"
  }
}


task ShapeIt4 {
  input {
    File merged_vcf
    File merged_vcf_index
    String project_id
    String contig
    File genetic_map
    Float window = 2.5
    Int num_threads = 2
    String mcmc_iterations = "5b,1p,1b,1p,1b,1p,5m"
    Int pbwt_depth = 4
    # TODO - how to handle refence as an option
    ReferenceSequence? reference

    String singularity_image = "shapeit4_4.1.3.sif"
    Int num_cpu = 4
    Int memory = 15000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    #TODO - handle reference...
    touch ~{merged_vcf_index}
    shapeit4 \
        --input ~{merged_vcf} \
        ~{"--map " + genetic_map} \
        --region ~{contig} \
        ~{"--window " + window} \
        ~{"--thread " + num_threads} \
        ~{"--mcmc-iterations " + mcmc_iterations} \
        ~{"--pbwt-depth " + pbwt_depth} \
        --sequencing \
        --use-PS 0.0001 \
        --log phased.log \
        --output ~{project_id}_phased.vcf.gz
  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
      File phased_vcf =  "~{project_id}_phased.vcf.gz"
      File log_file =  "phased.log"
  }
}

task CohortVcfToZarr {
  input {
    File input_vcf
    String contig
    String output_zarr_file_name
    String output_log_file_name

    String singularity_image = "cohortvcftozarr.1.1.sif"
    Int num_cpu = 1
    Int memory = 7500
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    mkdir outputs
    python /tools/cohort_vcf_to_zarr.py \
    --input ~{input_vcf} \
    --output ~{output_zarr_file_name} \
    --contig ~{contig} \
    --field variants/POS \
    --field variants/REF:S1 \
    --field variants/ALT:S1 \
    --field variants/AC \
    --field variants/AF \
    --field variants/CM \
    --field calldata/GT \
    --log ~{output_log_file_name} \
    --zip
  }
  runtime {
    singularity: singularity_image
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