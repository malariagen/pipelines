version 1.0

import "../../structs/gcp/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task MergeVcfs {
  input {
    Array[File] phased_sample_vcfs
    Array[File] phased_sample_vcf_indices
    String project_id

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/bcftools:1.11"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
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
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
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
    String region
    File genetic_map
    Float window = 2.5
    Int num_threads = num_cpu
    String mcmc_iterations = "5b,1p,1b,1p,1b,1p,5m"
    Int pbwt_depth = 4
    # TODO - how to handle refence as an option
    ReferenceSequence? reference

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/shapeit4:4.1.3"
    # Compute Engine always stops preemptible instances after they run for 24 hours
    Int preemptible_tries = 0
    Int num_cpu = 4
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(merged_vcf, "GiB") + size(merged_vcf_index, "GiB") + size(genetic_map, "GiB")) * 2 + 20
  String output_prefix = sub(region, ":", "_")
  String output_filename = output_prefix + "_" +  project_id + "_phased.vcf.gz"

  command {
    #TODO - handle reference...
    touch ~{merged_vcf_index}
    shapeit4 \
        --input ~{merged_vcf} \
        ~{"--map " + genetic_map} \
        --region ~{region} \
        ~{"--window " + window} \
        ~{"--thread " + num_threads} \
        ~{"--mcmc-iterations " + mcmc_iterations} \
        ~{"--pbwt-depth " + pbwt_depth} \
        --sequencing \
        --use-PS 0.0001 \
        --log phased.log \
        --output ~{output_filename}
  }

  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "15 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
      File region_phased_vcf =  "~{output_filename}"
      File log_file =  "phased.log"
  }
}



task LigateRegions {
  input {
    Array[File] region_phased_vcfs
    Array[File] region_phased_vcfs_indices
    File interval_list
    String project_id

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/bcftoolspython:1.11"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
  }

  Int disk_size = ceil(size(region_phased_vcfs, "GiB")) * 2 + 20

  command {
    set -e pipefail

    python3 <<CODE
    from pathlib import Path

    index_files = [ "~{sep='", "' region_phased_vcfs_indices}" ]
    for f in index_files:
        Path(f).touch()

    input_files = [ "~{sep='", "' region_phased_vcfs}" ]
    with open("~{interval_list}") as f:
        intervals = [i.replace(":", "_").strip("\n") for i in f]
    ordered_files = []

    for i in intervals:
        for f in input_files:
            if i in f:
                ordered_files.append(f)

    with open("phased_vcf_list.txt", "w") as f:
        f.write("\n".join(ordered_files))
    CODE

    bcftools concat \
        --file-list "phased_vcf_list.txt" \
        --ligate \
        --output ~{project_id}_phased.vcf

    bgzip -c ~{project_id}_phased.vcf > ~{project_id}_phased.vcf.gz
    tabix ~{project_id}_phased.vcf.gz
  }

  runtime {
    docker: docker_tag
    preemptible: preemptible_tries
    cpu: num_cpu
    memory: "15 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
      File phased_vcf =  "~{project_id}_phased.vcf.gz"
      File phased_vcf_index = "~{project_id}_phased.vcf.gz.tbi"
  }
}


task CohortVcfToZarr {
  input {
    File input_vcf
    String contig
    String output_zarr_file_name
    String output_log_file_name

    String docker_tag = "us.gcr.io/broad-gotc-prod/malariagen/cohortvcftozarr:1.1"
    Int preemptible_tries = runTimeSettings.preemptible_tries
    Int num_cpu = 1
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_vcf, "GiB")) * 4) + 20

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