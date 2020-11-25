version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"

task ImportIRODS {
  input {
    String irods_path
    String sample_id

    String singularity_image = runTimeSettings.irods_singularity_image
    Int num_cpu = 1
    Int memory = 1000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  String lanelet_basename = basename(irods_path)

  command {

    set -e
    set -o pipefail


    # Test that we can download from IRODS
    # -K verifies checksum
    # -f force overwite
    iget -K -f ${irods_path} ./${sample_id}.${lanelet_basename}
  }

  runtime {
    singularity: singularity_image
    memory: memory
    cpu: num_cpu
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
      File output_file = "${sample_id}.${lanelet_basename}"
  }

}


# Given a tsv sample manifest where each row contains
# file locations of each lanelet of sequencing,
# separate each sample into separate tsv manifests with no header.
# The per-sample manifests can then be used to process each sample in parallel tasks,
# where each task in turn scatters the per-sample manifests by lanelet.

task BatchSplitUpInputFile {

  input {
    File batch_sample_manifest_file

    String singularity_image = runTimeSettings.binder_singularity_image
    Int num_cpu = 1
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command <<<

  python <<CODE
  import pandas as pd

  sample_mf = pd.read_csv("~{batch_sample_manifest_file}", sep="\t")
  sample_mf = sample_mf.sort_values(["sample_id", "irods_path"])

  assert 'sample_id' in sample_mf.columns, "sample_id is not in header"
  assert 'irods_path' in sample_mf.columns, "irods_path is not in header"

  # Split the manifest into separate files
  for sample_id in sample_mf["sample_id"].unique():
      sample_mf[sample_mf["sample_id"] == sample_id].to_csv("{}_manifest.tsv".format(sample_id),
                                                            header=False,
                                                            sep="\t",
                                                            index=False)
  CODE
  >>>

  runtime {
    singularity: singularity_image 
    memory: memory
    cpu: num_cpu
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
      Array[File] per_sample_manifest_files = glob("*_manifest.tsv")
  }

}
