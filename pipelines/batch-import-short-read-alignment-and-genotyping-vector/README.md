# BatchImportShortReadAlignmentAndGenotyping (vector) pipeline

This folder contains the pipeline wdl and workflow inputs for the combined Vector Short Read Alignment Pipeline, SNP Genotyping pipeline as specified in https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md and https://github.com/malariagen/pipelines/blob/991b0328ea8027bc6b1137f893a5340e27c8e87c/docs/specs/snp-genotyping-vector.md


There are farm5-specific directories which contain versions specific to farm5 (Sanger LSF) backed instances of cromwell.  These pipelines import multiple lanelet bams/crams from Sanger IRODS file system for multiple samples in parallel.

The BatchImportShortReadAlignmentAndGenotyping pipeline takes as input:
- batch_sample_manifest_file: A tab-delimited file with columns for sample_id, irods_path.  Detailed below
- known_indels_vcf: An optional file of known indels for the tasks RealignerTargetCreator and IndelRealigner
- alleles_vcf: The alleles VCF for use by the task UnifiedGenotyper
- alleles_vcf_index: The index file of the alleles VCF
- reference: A wdl structure containing Vector reference files
- runTimeSettings: A wdl structure containing run time settings.

batch_sample_manifest_file: The input file is a tab-delimited file.  There is a header, and after that, each line represents a 'lanelet'.  The file can contain lanelets for multiple samples.  The pipeline will separate each sample into separate manifest files to be processed in parallel.


The mandatory fields are:
- sample_id: The Id/name of the sample.
- irods_path: The irods path for the lanelet data.  This can be a cram or a bam.  The basename will be used as the read group ID in the output aligned bam.

Additional fields may be specified, but they will not be used directly in the pipeline.

Example File:
The first couple of lines for an example input file are specified here:

| sample_id 	| irods_path                |
|-----------	|---------------------------|
| AN0131-C  	| /seq/9812/9812_4#48.bam   	|
| AN0131-C  	| /seq/10209/10209_3#48.bam   |
| AN0131-C  	| /seq/10209/10209_4#48.bam   |
| AB0252-C  	| /seq/9953/9953_2#68.bam     |
| AB0252-C  	| /seq/10061/10061_1#68.bam 	|
