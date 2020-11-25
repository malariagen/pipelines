version 1.0

import "../../structs/farm5/RunTimeSettings.wdl"
import "../../structs/ReferenceSequence.wdl"

task SplitUpInputFile {
  input {
    File input_file
    String sample_id

    String singularity_image = runTimeSettings.lftp_singularity_image
    Int num_cpu = 1
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  String lanelet_file_prefix = "lanelet_"

  command {
    mkdir lanelet_temp
    cd lanelet_temp

    # Verify that the input_file begins with a header
    # should look like: `sample_id	run_ena	irods_path	bam_path	cram_path	read1_path	read2_path`
    head -1 ~{input_file} | grep '^sample_id\trun_ena'
    exitCode=$?
    if [ $exitCode != 0 ]; then
      echo "Input file ~{input_file} appears malformed"
      exit $exitCode
    fi

    # splits list of mappings into single files.  One line each.
    grep '^~{sample_id}\t' ~{input_file} | split -l 1 - ~{lanelet_file_prefix}
  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
      Array[File] lanelet_files = glob("lanelet_temp/~{lanelet_file_prefix}*")
  }
}

task Ftp {
  input {
    String input_string
    String output_filename = basename(input_string)

    String singularity_image = runTimeSettings.lftp_singularity_image
    Int num_cpu = 1
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    echo get1 ~{input_string} -o ~{output_filename} > script_file
    lftp -f script_file

  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }

  output {
    File output_file = output_filename
  }
}

task CramToBam {
  input {
    File input_file
    String output_filename

    String singularity_image = runTimeSettings.samtools_singularity_image
    Int num_cpu = 2
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx~{memory}m -jar /bin/picard.jar \
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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx~{memory}m -jar /bin/picard.jar \
      SamToFastq \
      INPUT=~{input_file} \
      FASTQ=~{output_fastq1_filename} \
      SECOND_END_FASTQ=~{output_fastq2_filename} \
      NON_PF=true
  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.bwa_singularity_image
    Int num_cpu = 4
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
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
    /bwa/bwa mem -M -K 100000000 -t 4 -T 0 -R '@RG\tID:~{read_group_id}\tSM:~{sample_id}\tCN:SC\tPL:ILLUMINA' ~{reference.ref_fasta} ~{fastq1} ~{fastq2} > ~{output_sam_basename}.sam
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_sam = "~{output_sam_basename}.sam"
  }
}

task ReadAlignmentPostProcessing {
  input {
    File input_sam
    String output_bam_basename

    String singularity_image = runTimeSettings.samtools_singularity_image
    Int num_cpu = 2
    Int memory = 3000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    /bin/samtools view -bu ~{input_sam} |
    /bin/samtools sort -n - |
    /bin/samtools fixmate - - |
    /bin/samtools sort - > ~{output_bam_basename}.bam

  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task SetNmMdAndUqTags {
  input {
    File input_bam
    String output_bam_basename

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx~{memory}m -jar /bin/picard.jar \
      SetNmMdAndUqTags \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference.ref_fasta} \
      IS_BISULFITE_SEQUENCE=false
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task MergeSamFiles {
  input {
    Array[File] input_files
    String output_filename

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx~{memory}m -jar /bin/picard.jar \
      MergeSamFiles \
      INPUT=~{sep=' INPUT=' input_files} \
      OUTPUT=~{output_filename}
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_file = output_filename
  }
}

task MarkDuplicates {
  input {
    File input_bam
    String output_filename

    String singularity_image = runTimeSettings.biobambam_singularity_image
    Int num_cpu = 1
    Int memory = 2000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
    /usr/local/bin/bammarkduplicates I=~{input_bam} O=~{output_filename} index=1
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.gatk_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx~{memory}m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T RealignerTargetCreator \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -o ~{output_interval_list_filename}
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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
    String output_bam_filename

    String singularity_image = runTimeSettings.gatk_singularity_image
    Int num_cpu = 1
    Int memory = 8000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx~{memory}m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T IndelRealigner \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          ~{"-known " + known_indels_vcf} \
          -targetIntervals ~{interval_list_file} \
          -o ~{output_bam_filename}
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File output_bam = output_bam_filename
  }
}

task FixMateInformation {
  input {
    File input_file
    String output_bam_basename

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 8000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {
   set -e
   set -o pipefail

    java -Xmx~{memory}m -jar /bin/picard.jar \
      FixMateInformation \
      INPUT=~{input_file} \
      OUTPUT=~{output_bam_basename}.bam \
      MAX_RECORDS_IN_RAM=300000 \
      CREATE_INDEX=true

      # FixMateInformation creates the bam index as foo.bai, we move it to foo.bam.bai
      mv ~{output_bam_basename}.bai ~{output_bam_basename}.bam.bai
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.picard_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -Xmx~{memory}m -jar /bin/picard.jar \
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
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File report_file = report_filename
  }
}

task SamtoolsStats {
  input {
    File input_file
    String report_filename

    String singularity_image = runTimeSettings.samtools_singularity_image
    Int num_cpu = 1
    Int memory = 2000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    /bin/samtools stats -r ~{reference.ref_fasta} ~{input_file} > ~{report_filename}

  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.samtools_singularity_image
    Int num_cpu = 1
    Int memory = 2000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    /bin/samtools idxstats ~{input_bam} > ~{report_filename}

  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File report_file = report_filename
  }
}

task SamtoolsFlagStat {
  input {
    File input_bam
    String report_filename

    String singularity_image = runTimeSettings.samtools_singularity_image
    Int num_cpu = 1
    Int memory = 2000
    String? lsf_group
    String? lsf_queue
    RunTimeSettings runTimeSettings
  }

  command {

    set -e
    set -o pipefail

    /bin/samtools flagstat ~{input_bam} > ~{report_filename}

  }

  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
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

    String singularity_image = runTimeSettings.gatk_singularity_image
    Int num_cpu = 1
    Int memory = 4000
    String? lsf_group
    String? lsf_queue
    ReferenceSequence reference
    RunTimeSettings runTimeSettings
  }

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx~{memory}m \
          -jar /usr/GenomeAnalysisTK.jar \
          -T CallableLoci \
          -I ~{input_bam} \
          -R ~{reference.ref_fasta} \
          --summary ~{summary_filename} \
          --minDepth 5
  }
  runtime {
    singularity: singularity_image
    cpu: num_cpu
    memory: memory
    lsf_group: select_first([runTimeSettings.lsf_group, lsf_group, "pathdev"])
    lsf_queue: select_first([runTimeSettings.lsf_queue, lsf_queue, "normal"])
  }
  output {
    File summary_file = summary_filename
  }
}
