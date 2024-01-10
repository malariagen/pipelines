version 1.0

task ConsolidateHMMOutput {
    meta {
        description: "This task takes the output from the HMM subpipeline and combines it into a single tarball."
    }
    parameter_meta {
        hmm_tarballs: "The output files from the HMM sub-pipeline. This is an array of tarballs, one for each sample."
        output_dir: "The output directory for the tarball."
    }
    input {
        Array[File] hmm_tarballs
        String output_dir
        # runtime values
        String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
        Int ram = if (ceil(size(hmm_tarballs, "MiB") * 5) + 8000) < 120000 then (ceil(size(hmm_tarballs, "MiB") * 5) + 8000) else 120000
        Int cpu = 16
        Int disk =  ceil(size(hmm_tarballs, "GiB") * 10) + 50
        Int preemptible = 3
        String runtime_zones
    }
    command <<<
        set -x
        echo "Current directory: "
        pwd
        #For each file in hmm_tarballs, extract the tarball and move the contents to the output directory
        for tarball in ~{sep=' ' hmm_tarballs} ; do
        echo "Extracting tarball: " $tarball
        tar --backup=numbered -zxvf $tarball
        ls -lht
        done
        ls -lht
        tar -zcvf ~{output_dir}.tar.gz ~{output_dir}
    >>>
    runtime {
        docker: docker
        memory: "${ram} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
        zones: runtime_zones
    }
    output {
        File consolidated_gz = "~{output_dir}.tar.gz"
    }
}


task CreateSampleManifest {
    meta {
        description: "This task creates a new line separated sample manifest, where each sample id is a separate line."
    }
    parameter_meta {
        sample_ids: "An Array of sample ids for all samples bing run together"
    }

    input {
        Array[String] sample_ids
        # runtime values
        String docker = "ubuntu:18.04"
        String ram = "3000 MiB"
        Int cpu = 1
        Int disk = 10
        Int preemptible = 3
        String runtime_zones
    }

    command <<<
       # I don't think we actually need to do anythign here
       echo "writing sample ids to output file"
    >>>

    runtime {
        docker: docker
        memory: ram
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
        zones: runtime_zones
    }
    output {
        File sample_manifest = write_lines(sample_ids)
    }
}