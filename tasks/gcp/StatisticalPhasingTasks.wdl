version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {
    Array[File] sample_phased_vcfs
    String project_id
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }
  command {

  }
  runtime {
    docker: runTimeSettings.gatk_docker
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
    String project_id
    File genetic_map
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    shapeit4 \
        --input ~{merged_vcf} \
        --map ~{genetic_map} \
        --region REGION \
        --window WINDOW \
        --thread N_THREADS \
        --mcmc-iterations MCMC_ITERATIONS \
        --pbwt-depth PBWT_DEPTH \
        --window WINDOW \
        --sequencing \
        --reference REFERENCE.vcf \
        --use-PS 0.0001 \
        --log phased.log \
        --output ~{project_id}_phased.vcf.gz
  }

  runtime {
    docker: runTimeSettings.shapeit4
    preemptible: runTimeSettings.preemptible_tries
    cpu: "1"
    memory: "3.75 GiB"
  }

  output {
      File phased_vcf =  "~{project_id}_phased.vcf.gz"
      Array[File] phasing_logs = glob("lanelet_temp/~{lanelet_file_prefix}*")
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