version 1.0

## Copyright Wellcome Sanger Institute, Oxford University, and the Broad Institute 2020
##
## This WDL pipeline implements the Short Read Alignment Pipeline as described in
## https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md
## This is an initial proof of concept implementation.  It is designed to ONLY work on one sample, with ONLY its
## list of input fastqs.  It is currently implemented to run using Cromwell with a google cloud platform backend.
##

workflow ShortReadAlignment {
  String pipeline_version = "0.0.1"

  input {
    String sample_id
    Array[String] read_group_ids
    Array[String] libraries
    Array[String] studies
    Array[File] input_fastq1s
    Array[File] input_fastq2s
    String output_basename

    File? known_indels_vcf

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  scatter (idx in range(length(input_fastq1s))) {
    call ReadAlignment {
      input:
        sample_id = sample_id,
        read_group_id = read_group_ids[idx],
        library = libraries[idx],
        study = studies[idx],
        fastq1 = input_fastq1s[idx],
        fastq2 = input_fastq2s[idx],
        output_sam_basename = output_basename,
        reference = reference,
        runTimeSettings = runTimeSettings
    }

    call ReadAlignmentPostProcessing {
      input:
        input_sam = ReadAlignment.output_sam,
        output_bam_basename = output_basename,
        runTimeSettings = runTimeSettings
    }

    call SetNmMdAndUqTags {
      input:
        input_bam = ReadAlignmentPostProcessing.output_bam,
        output_bam_basename = output_basename,
        reference = reference,
        runTimeSettings = runTimeSettings
    }
  }

  call MergeSamFiles {
    input:
      input_files = SetNmMdAndUqTags.output_bam,
      output_filename = output_basename + ".bam",
      runTimeSettings = runTimeSettings
  }

  call MarkDuplicates {
    input:
      input_bam = MergeSamFiles.output_file,
      output_filename = output_basename + ".bam",
      runTimeSettings = runTimeSettings
  }

  call RealignerTargetCreator {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      output_interval_list_filename = output_basename + ".intervals",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call IndelRealigner {
    input:
      input_bam = MarkDuplicates.output_file,
      input_bam_index = MarkDuplicates.output_index_file,
      interval_list_file = RealignerTargetCreator.output_interval_list_file,
      output_filename = output_basename + ".bam",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call FixMateInformation {
    input:
      input_file = IndelRealigner.output_file,
      output_bam_basename = output_basename,
      runTimeSettings = runTimeSettings
  }

  call ValidateSamFile {
    input:
      input_file = FixMateInformation.output_file,
      report_filename = output_basename + ".validation_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call SamtoolsStats {
    input:
      input_file = FixMateInformation.output_file,
      report_filename = output_basename + ".samtools_stats_report.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  call GatkCallableLoci {
    input:
      input_bam = FixMateInformation.output_file,
      input_bam_index = FixMateInformation.output_index_file,
      summary_filename = output_basename + ".gatk_callable_loci.summary.txt",
      reference = reference,
      runTimeSettings = runTimeSettings
  }

  output {
    Array[File] output_sams = ReadAlignment.output_sam
    Array[File] output_bams = ReadAlignmentPostProcessing.output_bam
    File indel_realigner_interval_list_file = RealignerTargetCreator.output_interval_list_file
    File output_bam = IndelRealigner.output_file
    File validate_samfile_report_file = ValidateSamFile.report_file
    File samtools_stats_report_file = SamtoolsStats.report_file
    File callable_loci_summary_file = GatkCallableLoci.summary_file
  }
}

task ReadAlignment {
  input {
    String read_group_id
    String library
    String sample_id
    String study
    File fastq1
    File fastq2
    String output_sam_basename

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Float fastq_size = size(fastq1, "GiB") + size(fastq2, "GiB")
  Float ref_size = size(reference.ref_fasta, "GiB") + size(reference.ref_fasta_index, "GiB") + size(reference.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference.ref_amb, "GiB") + size(reference.ref_ann, "GiB") + size(reference.ref_bwt, "GiB") + size(reference.ref_pac, "GiB") + size(reference.ref_sa, "GiB")
  Float disk_multiplier = 2.5
  Int disk_size = ceil(fastq_size + bwa_ref_size + (disk_multiplier * fastq_size) + 20)

  command <<<
    set -o pipefail
    set -e

    # suggested content for the read group tag: (from Thuy):
    # Read group identifier [ID]: full platform unit ID as the read group identifier (including the flowcell, run, lane)
    # library [LB]: obtained from raw sequenced bam, but we can make this up for testing
    # sample [SM]: obtained from raw sequenced bam, but we can also obtain this from the sample manifest or fastq filename
    # sequencing centre [CN]: obtained from raw sequenced bam, but we can make this up for testing
    # platform [PL]: obtained from raw sequenced bam, but we can make this up for testing
    # study [DS]: obtained from raw sequenced bam, but we can make this up for testing
    # @rg ID:130508_HS22_09812_A_D1U5TACXX_4#48 LB:7206533 SM:AN0131-C CN:SC PL:ILLUMINA DS:1087-AN-HAPMAP-DONNELLY
    /bwa/bwa mem -M -T 0 -R '@RG\tID:~{read_group_id}\tLB:~{library}\tSM:~{sample_id}\tCN:SC\tPL:ILLUMINA\tDS:~{study}' ~{reference.ref_fasta} ~{fastq1} ~{fastq2} > ~{output_sam_basename}.sam
  >>>
  runtime {
    docker: runTimeSettings.bwa_docker
    preemptible: runTimeSettings.preemptible_tries
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_sam = "~{output_sam_basename}.sam"
  }
}

task ReadAlignmentPostProcessing {
  input {
    File input_sam
    String output_bam_basename
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_sam, "GiB")) * 3) + 20

  command <<<

    set -e
    set -o pipefail

    /bin/samtools view -bu ~{input_sam} |
    /bin/samtools sort -n - |
    /bin/samtools fixmate - - |
    /bin/samtools sort - > ~{output_bam_basename}.bam

  >>>

  runtime {
    docker: runTimeSettings.samtools_docker
    preemptible: runTimeSettings.preemptible_tries
    memory: "14 GiB"
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task SetNmMdAndUqTags {
  input {
    File input_bam
    String output_bam_basename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 3) + 20

  command {
    java -Xms3500m -jar /bin/picard.jar \
      SetNmMdAndUqTags \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      IS_BISULFITE_SEQUENCE=false
  }
  runtime {
    docker: runTimeSettings.picard_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task MergeSamFiles {
  input {
    Array[File] input_files
    String output_filename
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_files, "GiB")) * 3) + 20

  command {
    java -Xms3500m -jar /bin/picard.jar \
      MergeSamFiles \
      INPUT=~{sep=' INPUT=' input_files} \
      OUTPUT=~{output_filename}
  }
  runtime {
    docker: runTimeSettings.picard_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
  }
  output {
    File output_file = output_filename
  }
}

task MarkDuplicates {
  input {
    File input_bam
    String output_filename
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 3) + 20

  command {
    /usr/local/bin/bammarkduplicates I=~{input_bam} O=~{output_filename} index=1
  }
  runtime {
    docker: runTimeSettings.biobambam_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
  }
  output {
    File output_file = output_filename
    File output_index_file = "~{output_filename}.bai"
  }
}

task RealignerTargetCreator {
  input {
    File input_bam
    File input_bam_index
    File? known_indels_vcf
    String output_interval_list_filename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 2) + 20

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T RealignerTargetCreator \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -o ~{output_interval_list_filename}
  }
  runtime {
    docker: runTimeSettings.gatk_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
  }
  output {
    File output_interval_list_file = output_interval_list_filename
  }
}

task IndelRealigner {
  input {
    File input_bam
    File input_bam_index
    File? known_indels_vcf
    File interval_list_file
    String output_filename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 2) + 20

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms7000m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T IndelRealigner \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -targetIntervals ~{interval_list_file} \
          -o ~{output_filename}
  }
  runtime {
    docker: runTimeSettings.gatk_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "7.5 GiB"
  }
  output {
    File output_file = output_filename
  }
}

task FixMateInformation {
  input {
    File input_file
    String output_bam_basename
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_file, "GiB")) * 3) + 20

  command {
    java -Xms7000m -jar /bin/picard.jar \
      FixMateInformation \
      INPUT=~{input_file} \
      OUTPUT=~{output_bam_basename}.bam \
      MAX_RECORDS_IN_RAM=300000 \
      CREATE_INDEX=true
  }
  runtime {
    docker: runTimeSettings.picard_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "7.5 GiB"
  }
  output {
    File output_file = "~{output_bam_basename}.bam"
    File output_index_file = "~{output_bam_basename}.bai"
  }
}

task ValidateSamFile {
  input {
    File input_file
    File? input_file_index
    String report_filename
    Int? max_output
    Array[String]? ignore
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Float ref_size = size(reference.ref_fasta, "GiB") + size(reference.ref_fasta_index, "GiB") + size(reference.ref_dict, "GiB")
  Int disk_size = ceil(size(input_file, "GiB") + ref_size) + 20

  command {
    java -Xms3500m -jar /bin/picard.jar \
      ValidateSamFile \
      INPUT=~{input_file} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      IS_BISULFITE_SEQUENCED=false \
      MODE=VERBOSE
  }
  runtime {
    docker: runTimeSettings.picard_docker
    preemptible: runTimeSettings.preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report_file = report_filename
  }
}

task SamtoolsStats {
  input {
    File input_file
    String report_filename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Float ref_size = size(reference.ref_fasta, "GiB") + size(reference.ref_fasta_index, "GiB") + size(reference.ref_dict, "GiB")
  Int disk_size = ceil(size(input_file, "GiB") + ref_size) + 20

  command <<<

    set -e
    set -o pipefail

    /bin/samtools stats -r ~{reference.ref_fasta} ~{input_file} > ~{report_filename}

  >>>

  runtime {
    docker: runTimeSettings.samtools_docker
    preemptible: runTimeSettings.preemptible_tries
    memory: "7.5 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report_file = report_filename
  }
}

task GatkCallableLoci {
  input {
    File input_bam
    File input_bam_index
    String summary_filename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  Int disk_size = (ceil(size(input_bam, "GiB")) * 2) + 20

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T CallableLoci \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          --summary ~{summary_filename} \
          --minDepth 5
  }
  runtime {
    docker: runTimeSettings.gatk_docker
    preemptible: runTimeSettings.preemptible_tries
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
  }
  output {
    File summary_file = summary_filename
  }
}


struct ReferenceSequence {
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File? ref_alt
  File ref_sa
  File ref_amb
  File ref_bwt
  File ref_ann
  File ref_pac
}

struct RunTimeSettings {
  Int? preemptible_tries
  String gatk_docker
  String picard_docker
  String bwa_docker
  String biobambam_docker
  String samtools_docker
}
