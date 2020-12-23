version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {
    Array[File] phased_sample_vcfs
    Array[File] phased_sample_vcf_indices
    String project_id
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(phased_sample_vcfs, "GiB") + size(phased_sample_vcf_indices, "GiB")) * 2 + 20

  #TODO - touch all the index files to avoid the annoying warning
  command {
    bcftools merge \
      -o ~{project_id}_merged.vcf \
      ~{sep=' ' phased_sample_vcfs}
  }
  runtime {
    docker: runTimeSettings.bcftools_docker
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
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

  Int disk_size = ceil(size(merged_vcf, "GiB") + size(merged_vcf_index, "GiB") + size(genetic_map, "GiB")) * 2 + 20

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
    cpu: "4"
    memory: "15 GiB"
    disks: "local-disk " + disk_size + " HDD"
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