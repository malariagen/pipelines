version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {
    Array[File] sample_phased_vcfs
    Array[File] sample_phased_vcf_indices
    String project_id
    RunTimeSettings runTimeSettings
  }
  command {
    bcftools merge \
      -o ~{project_id}_merged.vcf \
      ~{sep=' ' sample_phased_vcfs}
  }
  runtime {
    docker: runTimeSettings.bcftools_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
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
    docker: runTimeSettings.shapeit4_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "2"
    memory: "7.5 GiB"
  }

  output {
      File phased_vcf =  "~{project_id}_phased.vcf.gz"
      File log_file =  "phased.log"
  }
}

task LigateRegions {
  input {

  }
  command {

  }
  runtime {

  }
  output {

  }
}

# task VcfToZarr {
# TODO: can we reuse the vcf to zarr from the genotyping pipeline?
# }