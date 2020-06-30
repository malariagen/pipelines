version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task CramToBam {
  input {
    File input_file
    String output_filename
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    samtools view -h -T ~{reference.ref_fasta} ~{input_file} |
    samtools view -b -o ~{output_filename} -
    samtools index -b ~{output_filename}

  }

  runtime {
    singularity: runTimeSettings.samtools_singularity_image
    memory: 3000
    cpu: "2"
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }

  output {
    File output_file = output_filename
    File output_file_index = "~{output_filename}.bai"
  }
}

task RevertSam {
  input {
    File input_file
    String output_filename
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx3500m -jar /bin/picard.jar \
      RevertSam \
      INPUT=~{input_file} \
      OUTPUT=~{output_filename} \
      VALIDATION_STRINGENCY=LENIENT \
      ATTRIBUTE_TO_CLEAR=FT \
      ATTRIBUTE_TO_CLEAR=CO \
      ATTRIBUTE_TO_CLEAR=PA \
      ATTRIBUTE_TO_CLEAR=OA \
      ATTRIBUTE_TO_CLEAR=XA
  }

  runtime {
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }

  output {
    File output_file = output_filename
  }
}

task SamToFastq {
  input {
    File input_file
    String output_fastq1_filename
    String output_fastq2_filename
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx3500m -jar /bin/picard.jar \
      SamToFastq \
      INPUT=~{input_file} \
      FASTQ=~{output_fastq1_filename} \
      SECOND_END_FASTQ=~{output_fastq2_filename} \
      NON_PF=true
  }

  runtime {
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }

  output {
    File output_fastq1 = output_fastq1_filename
    File output_fastq2 = output_fastq2_filename
  }
}

task ReadAlignment {
  input {
    String read_group_id
    String sample_id
    File fastq1
    File fastq2
    String output_sam_basename

    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

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
    /bwa/bwa mem -M -T 0 -R '@RG\tID:~{read_group_id}\tSM:~{sample_id}\tCN:SC\tPL:ILLUMINA' ~{reference.ref_fasta} ~{fastq1} ~{fastq2} > ~{output_sam_basename}.sam
  >>>
  runtime {
      singularity: runTimeSettings.bwa_singularity_image
      memory: 1500
      cpu: "1"
      lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
      lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command <<<

    set -e
    set -o pipefail

    /bin/samtools view -bu ~{input_sam} |
    /bin/samtools sort -n - |
    /bin/samtools fixmate - - |
    /bin/samtools sort - > ~{output_bam_basename}.bam

  >>>

  runtime {
    singularity: runTimeSettings.samtools_singularity_image
    memory: 3000
    cpu: "2"
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    java -Xmx3500m -jar /bin/picard.jar \
      SetNmMdAndUqTags \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      IS_BISULFITE_SEQUENCE=false
  }
  runtime {
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    java -Xmx3500m -jar /bin/picard.jar \
      MergeSamFiles \
      INPUT=~{sep=' INPUT=' input_files} \
      OUTPUT=~{output_filename}
  }
  runtime {
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    /usr/local/bin/bammarkduplicates I=~{input_bam} O=~{output_filename} index=1
  }
  runtime {
    singularity: runTimeSettings.biobambam_singularity_image
    memory: 1000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T RealignerTargetCreator \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -o ~{output_interval_list_filename}
  }
  runtime {
    singularity: runTimeSettings.gatk_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx7500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T IndelRealigner \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -targetIntervals ~{interval_list_file} \
          -o ~{output_filename}
  }
  runtime {
    singularity: runTimeSettings.gatk_singularity_image
    memory: 8000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
   set -e
   set -o pipefail

    java -Xmx7000m -jar /bin/picard.jar \
      FixMateInformation \
      INPUT=~{input_file} \
      OUTPUT=~{output_bam_basename}.bam \
      MAX_RECORDS_IN_RAM=300000 \
      CREATE_INDEX=true

      # FixMateInformation creates the bam index as foo.bai, we move it to foo.bam.bai
      mv ~{output_bam_basename}.bai ~{output_bam_basename}.bam.bai
  }
  runtime {
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bam.bai"
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

  command {
    java -Xmx3500m -jar /bin/picard.jar \
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
    singularity: runTimeSettings.picard_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command <<<

    set -e
    set -o pipefail

    /bin/samtools stats -r ~{reference.ref_fasta} ~{input_file} > ~{report_filename}

  >>>

  runtime {
    singularity: runTimeSettings.samtools_singularity_image
    cpu: "1"
    memory: 1000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }
  output {
    File report_file = report_filename
  }
}

task SamtoolsIdxStats {
  input {
    File input_bam
    File input_bam_index
    String report_filename
    RunTimeSettings runTimeSettings
  }

  command <<<

    set -e
    set -o pipefail

    /bin/samtools idxstats ~{input_bam} > ~{report_filename}

  >>>

  runtime {
    singularity: runTimeSettings.samtools_singularity_image
    cpu: "1"
    memory: 1000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }
  output {
    File report_file = report_filename
  }
}

task SamtoolsFlagStat {
  input {
    File input_bam
    String report_filename
    RunTimeSettings runTimeSettings
  }

  command <<<

    set -e
    set -o pipefail

    /bin/samtools flagstat ~{input_bam} > ~{report_filename}

  >>>

  runtime {
    singularity: runTimeSettings.samtools_singularity_image
    cpu: "1"
    memory: 1000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
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

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3500m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T CallableLoci \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          --summary ~{summary_filename} \
          --minDepth 5
  }
  runtime {
    singularity: runTimeSettings.gatk_singularity_image
    memory: 4000
    lsf_group: select_first([runTimeSettings.lsf_group, "malaria-dk"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, "normal"])
  }
  output {
    File summary_file = summary_filename
  }
}