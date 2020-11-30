# ImportShortReadAlignment (vector) pipeline

This folder contains the pipeline wdl and workflow inputs for the Vector Short Read Alignment pipeline as specified in https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md

There are farm5-specific directories which contain versions specific to farm5 (Sanger LSF) backed instances of cromwell.  Imports lanelet bams/crams from Sanger IRODS file system.

The ShortReadAlignment pipeline takes as input:
- sample_id: The Id/name of the sample to analyze
- input_file: A tab-delimited file with columns for sample_id, and various possible sources of data for the pipeline.  Detailed below
- known_indels_vcf: An optional file of known indels for the tasks RealignerTargetCreator and IndelRealigner
- reference: A wdl structure containing Vector reference files
- runTimeSettings: A wdl structure containing run time settings.

input_file: The input file is a tab-delimited file.  There is a header, and after that, each line represents a 'lanelet'
The fields are:
- sample_id: The Id/name of the sample.  The pipeline will only parse lines out of the input_file where the sample_id matches the wdl input 'sample_id'
- run_ena: The ena for the lanelet
- irods_path: The irods path for the lanelet data.  This can be a cram or a bam


per_sample_manifest_file: The input file is a tab-delimited file.  There is a header, and after that, each line represents a 'lanelet'.  Only the lanelets for a single sample should be specified.
The fields are:
- sample_id: The Id/name of the sample.  The pipeline will only parse lines out of the input_file where the sample_id matches the wdl input 'sample_id'
- run_ena: The ena for the lanelet.  At the moment, this is only used to name the output of alignment for each lanelet.  
- irods_path: The irods path for the lanelet data.  This can be a cram or a bam.

Example File:
The first couple of lines for an example input file are specified here:

| sample_id 	| run_ena   	| irods_path                	|
|-----------	|-----------	|---------------------------	|
| AN0131-C  	| ERR317337 	| /seq/9812/9812_4#48.bam   	|
| AN0131-C  	| ERR340933 	| /seq/10209/10209_3#48.bam 	|
| AN0131-C  	| ERR340945 	| /seq/10209/10209_4#48.bam 	|
