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
        Array[String] dependency_string
        # runtime values
        String docker = "us.gcr.io/broad-gotc-prod/cnv:1.0.0-1679431881"
        String ram = "8000 MiB"
        Int cpu = 16
        # TODO: Make disk space dynamic based on input size
        Int disk = 70
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
        memory: ram
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
        zones: runtime_zones
    }
    output {
        File consolidated_gz = "~{output_dir}.tar.gz"
    }
}

task GetSpeciesManifests {
    meta {
        description: "This task takes the species_id_file and creates species specific manifest files."
    }
    parameter_meta {
        species_id_file: "The samples.species_aim.csv (may also be in tsv format) that includes species information for each sample."
    }

    input {
        File species_id_file
        # runtime values
        String docker = "ubuntu:18.04"
        String ram = "3000 MiB"
        Int cpu = 1
        Int disk = 10
        Int preemptible = 3
    }

    command <<<
        # Need to create species specific manifests
        # Not clear whether we will get the species_id_file in csv or tsv format, so we'll convert just to be safe
        sed -i -e "s/,/\t/g" ~{species_id_file}

        # We don't want to assume a column ordering since this is not gaurenteed, instead we will parse by column name
        ID_COL=$(awk -v RS='\t' ''/^sample_id$/{print NR; exit}' ~{species_id_file})
        SPECIES_COL=$(awk -v RS='\t' '/^aim_species_gambcolu_arabiensis$/{print NR; exit}' ~{species_id_file})

        cut -f $ID_COL,$SPECIES_COL  ~{species_id_file} | tail -n +2 > species_manifest.tsv
        cut -f 2 species_manifest.tsv  | sort -u > all_species.txt
        cat all_species.txt | while read species; do grep $species species_manifest.tsv | cut -f 1 > sample_manifest_$species.txt; done
    >>>

    runtime {
        docker: docker
        memory: ram
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }
    output {
        Array[File] species_specific_sample_manifests = "sample_manifest_*"
        Array[String] species = read_lines("all_species.txt")
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