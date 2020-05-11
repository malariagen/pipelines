# Short read alignment (vector) pipeline specification

* Version: 1.2.1
* Authors: Alistair Miles, Jim Stalker, George Grant

This document specifies a protocol for alignment of short sequence
reads, intended for use with mosquito specimens.


## Assumptions

* The pipeline will be launched to align sequence reads from one or
  more samples as a batch.

* Any demultiplexing of multiple samples within a single sequencing
  lane will have been resolved upstream of this pipeline. I.e., the
  pipeline will run from demultiplexed sequence reads. The term
  "lanelet" is used here to mean data from a single lane of sequencing
  from a single sample.

* Each sample may have data from multiple sequencing runs (lanes), and
  the number of lanelets may vary between samples. I.e., the pipeline
  cannot assume a fixed number of lanelets per sample.

* The complete set of input data will be known ahead of time, before
  the pipeline is launched for a given batch of samples. I.e., there
  is no need to consider any case where additional lanelets may arrive
  at a later time.

* Sequencing will be from Illumina technology, using paired-end
  sequencing. Mosquitoes will be sequenced individually, i.e., each
  mosquito sample will be used to create a separate DNA library, and
  libraries will not be pooled. Each individual mosquito will be
  sequenced to high coverage (target ~30X). Read length may vary
  between sequencing runs, i.e., the pipeline cannot assume the same
  read length across all lanelets.

* Trimming of adapter sequences, if necessary, will have been
  performed prior to running this pipeline.


## Inputs

* Reference sequence - The reference genome sequence to align to, in
  FASTA format.

* Sequence reads - For each lanelet, a pair of files containing
  sequence reads, with one file containing forward reads and one file
  containing reverse reads, in gzipped FASTQ format.

* Lanelet manifest - A tab delimited file with one header row followed
  by one data row per lanelet. Mandatory columns are as follows and
  are expected to occur as the first columns in the file. Other
  columns may be included (will be ignored).

    * `read1_path` - URL or local filesystem path to gzipped FASTQ
      file containing forward reads.

    * `read2_path` - URL or local filesystem path to gzipped FASTQ
      file containing reverse reads.

    * `sample_id` - Identifier for the sample that was sequenced.


## Outputs

* One analysis-ready BAM file for each sample.

* samtools stats output for each sample in TSV format.

* GATK CallableLoci output for each sample.


## Software

* bwa 0.7.15
* samtools 1.4.1 (htslib 1.4.1)
* picard 2.9.2
* biobambam 2.0.73
* GATK 3.7-0


## Pipeline


### Sub-pipeline: Lanelet alignment

For each lanelet, the following steps are run.


#### Step: Read alignment

Align sequence reads for a single lanelet to the reference sequence with command:

```
bwa mem -M -T 0 -R readgroup_string ref.fa read1.fastq read2.fastq > output.sam
```

Notes:

* -M marks split alignments as secondary (for Picard compatibility)

* -T 0 sets the threshold for alignment score output to zero


#### Step: Read alignment post-processing

Conversion of sam to bam via namesort, fixmates, and coordinate sort:

```
samtools view -bu output.sam 
| samtools sort -n -
| samtools fixmate - - 
| samtools sort - > output.bam
```

#### Step: Read alignment post-processing

Calculate the NM tags in the bam:

```
picard.jar SetNmMdAndUqTags input.bam
```


### Sub-pipeline: Sample alignment merging and improvement

For each sample, the following steps are run.


#### Merge to sample

Merge alignments from all lanes to build a BAM/CRAM for a single sample:

```
picard.jar MergeSamFiles lane1.bam lane2.bam lane3.bam 
```


#### Mark duplicate reads

Mark duplicate reads with command `bammarkduplicates` from biobambam.


#### Indel realignment

Realign reads around indels found within the alignments using the following commands:

```
GATK -T RealignerTargetCreator
```

Note that the RealignerTargetCreator tool has an option to input a
list of any known indels. However, we do not have any prior list of
known indels for our mosquito projects, and so we don't use this
option.

```
GATK -T IndelRealigner -targetIntervals output_of_RTCreator
```

Default options for LOD, model, etc.


#### Fix Mates

Correct mate pair information in the bam:

```
picard.jar FixMateInformation
```


#### Validation

Basic check for bam file validity, as interpreted by the Broad.

```
picard.jar ValidateSamFile
```

At this point, the analysis-ready BAM has been produced and validated.
Subsequent steps produce metrics from the BAM.


#### Samtools stats

Collects statistics from BAM files

```
samtools stats -r ref.fa input.bam
```


#### GATK callable loci

Collects statistics on callable, uncallable, poorly mapped, and other parts of the genome

```
GATK -T CallableLoci -R ref.fa -I input.bam --minDepth 5 
```


## Implementation notes

To achieve parallelisation, the lanelet alignment sub-pipeline can be run
in parallel for each lanelet. Some further parallelisation could also
be achieved within the lanelet pipeline by splitting the input files
into chunks, processing each chunk in parallel, then merging chunks
back together.

The sample alignment merging and improvement sub-pipeline can be run
in parallel for each sample. Some further parallelisation could also
be achieved within this sub-pipeline by splitting the input file into
chunks by genome region, processing each chunk in parallel, then
merging chunks back together.


## Change log

* Version 1.2.1 - Replace Samtool's calmd with Picard's SetNmMdAndUqTags for calculation of NM tags

* Version 1.1.2 - Small correction to the example samtools commands in the read alignment postprocessing step.

* Version 1.1.1 - Some editing of the spec to improve clarity,
  particularly around how lanelets and samples are processed.

* Version 1.1.0 - Version ported across from MalariaGEN vector-ops
  repo.
