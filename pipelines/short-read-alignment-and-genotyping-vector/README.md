# ShortReadAlignmentAndGenotyping (vector) pipeline

This folder contains the pipeline wdl and workflow inputs for the combined Vector Short Read Alignment Pipeline and SNP Genotyping pipeline as specified in https://github.com/malariagen/pipelines/blob/c7210d93628aaa31f26baa88a92e10322368b78e/docs/specs/short-read-alignment-vector.md and https://github.com/malariagen/pipelines/blob/991b0328ea8027bc6b1137f893a5340e27c8e87c/docs/specs/snp-genotyping-vector.md

There are gcp and farm5-specific directories which contain versions specific to gcp (Google Cloud Platform) and farm5 (Sanger LSF) backed instances of cromwell.  For the wdl (pipelines) the ONLY differences are in the import statements.  For the example inputs the paths in those are specific to the platform in use.

The ShortReadAlignmentAndGenotyping pipeline takes as input:
- sample_id: The Id/name of the sample to analyze
- output_base_name: The base filename that can be given to the outputs of the pipeline.  If not specified, the sample_id is used.
- input_file: A tab-delimited file with columns for sample_id, and various possible sources of data for the pipeline.  Detailed below
- known_indels_vcf: An optional file of known indels for the tasks RealignerTargetCreator and IndelRealigner
- alleles_vcf: The alleles VCF for use by the task UnifiedGenotyper
- alleles_vcf_index: The index file of the alleles VCF
- reference: A wdl structure containing Vector reference files
- runTimeSettings: A wdl structure containing run time settings.

input_file: The input file is a tab-delimited file.  There is a header, and after that, each line represents a 'lanelet'
The fields are:
- sample_id: The Id/name of the sample.  The pipeline will only parse lines out of the input_file where the sample_id matches the wdl input 'sample_id'
- run_ena: The ena for the lanelet
- irods_path: The irods path for the lanelet data.  This can be a cram or a bam
- bam_path: The bam path for the lanelet data.  If this is an ftp path (begins with ftp://) the pipeline will download the data from the ftp source
- cram_path: The cram path for the lanelet data.  If this is an ftp path (begins with ftp://) the pipeline will download the data from the ftp source
- read1_path: The path to the fastq1 data file for the lanelet.  If this is an ftp path (begins with ftp://) the pipeline will download the data from the ftp source
- read2_path: The path to the fastq2 data file for the lanelet.  If this is an ftp path (begins with ftp://) the pipeline will download the data from the ftp source

Note: If there is data specified in the cram_path column, that data will be used.  If there is none, then bam_path will be preferred, failing that the pipeline will use read1_path and read2_path for lanelet data.

Example File:
The first couple of lines for an example input file are specified here:

| sample_id 	| run_ena   	| irods_path                	| bam_path                                                           	| cram_path 	| read1_path                                                               	| read2_path                                                               	|
|-----------	|-----------	|---------------------------	|--------------------------------------------------------------------	|-----------	|--------------------------------------------------------------------------	|--------------------------------------------------------------------------	|
| AN0131-C  	| ERR317337 	| /seq/9812/9812_4#48.bam   	| ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR317/ERR317337/9812_4%2348.bam  	|           	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR317/ERR317337/ERR317337_1.fastq.gz 	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR317/ERR317337/ERR317337_2.fastq.gz 	|
| AN0131-C  	| ERR340933 	| /seq/10209/10209_3#48.bam 	| ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR340/ERR340933/10209_3%2348.bam 	|           	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR340/ERR340933/ERR340933_1.fastq.gz 	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR340/ERR340933/ERR340933_2.fastq.gz 	|
| AN0131-C  	| ERR340945 	| /seq/10209/10209_4#48.bam 	| ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR340/ERR340945/10209_4%2348.bam 	|           	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR340/ERR340945/ERR340945_1.fastq.gz 	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR340/ERR340945/ERR340945_2.fastq.gz 	|
| AB0252-C  	| ERR327108 	| /seq/9953/9953_2#68.bam   	| ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR327/ERR327108/9953_2%2368.bam  	|           	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR327/ERR327108/ERR327108_1.fastq.gz 	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR327/ERR327108/ERR327108_2.fastq.gz 	|
| AB0252-C  	| ERR338386 	| /seq/10061/10061_1#68.bam 	| ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR338/ERR338386/10061_1%2368.bam 	|           	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338386/ERR338386_1.fastq.gz 	| ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338386/ERR338386_2.fastq.gz 	|