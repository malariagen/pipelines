
# Dockerfile images

### amplicon-parasite-tools
Contains a conda environment with all tools needed for the amplicon SNP-calling (BCL-to-VCF, "Stage1, Steps 1-3") parasite pipeline
(except for those tools used to process manifest files, see [here](https://gitlab.com/malariagen-aspis/aspis-pipeline/-/blob/e8b35283ad70c41b4d0c9c9a3a660e31fff4431b/Docker/ManifestTools/Dockerfile)).
Environment file copied from [here](https://gitlab.com/malariagen-aspis/aspis-pipeline/-/blob/e8b35283ad70c41b4d0c9c9a3a660e31fff4431b/Stage1-Pipeline/conda/pipe-tools.txt).

### BWA
Supports the `build-arg` `version` that represents the tag name found in the repo
https://github.com/lh3/bwa. Default version v0.7.17.

### CohortVcfToZarr
Contains the cohort_vcf_to_zarr python script

### Samtools
Supports the `build-arg` `version` that represents the tag name found in the repo
https://github.com/samtools/samtools. Default version 1.9.

### BioBamBam
Comes in two different flavours v2.0.73 and v2.0.106. New versions are quick to build but can not be automated due to the release methods used for this program.

### Picard
Supports the `build-arg` `version` that represents the tag name found in the repo
https://github.com/broadinstitute/picard. Default version 2.22.3. Executed using `java -jar /bin/picard.jar`

### GATK
The Broad Institute provided Dockerfile for GATK can be found on [docker hub](https://hub.docker.com).
- For GATK 3.X use https://hub.docker.com/r/broadinstitute/gatk3
- For GATK 4.X use https://hub.docker.com/r/broadinstitute/gatk

The GATK 3.7-0 image (provided by the Broad Institute) has been been extended with htslib software for use in Sanger pipelines.

### Lftp
Contains the lftp too, a simple ftp client, used for downloading files from ena

### SampleVcfToZarr
Contains the sample_vcf_to_zarr python script

### SampleSelectVariants
Contains the sample_select_variants python script

### WhatsHap
Contains the latest version of WhatsHap https://whatshap.readthedocs.io/en/latest/

### Shapeit4
Contains the latest version of Shapeit4 https://odelaneau.github.io/shapeit4/

### Bcftools
Supports the `build-arg` `version` that represents the tag name found in the repo
https://github.com/bcftools/bcftools. Default version 1.11.

### Import
A farm5-specific docker image, providing dependencies for the iRODS import and batch processing workflows.

## Mosquito short read alignment, SNP Genotyping and Phasing pipelines

- bwa 0.7.15 - `docker build . --tag=bwa:0.7.15 --build-arg version=v0.7.15`
- samtools 1.4.1 (htslib 1.4.1) - `docker build . --tag=samtools:1.4.1 --build-arg version=1.4.1`
- picard 2.9.2 - ` docker build . --tag=picard:2.9.2 --build-arg version=2.9.2`
- biobambam 2.0.73 - `docker build . --tag=biobambam:2.0.73`
- GATK 3.7-0 - `docker build . --tag=sangerpathogens/gatk3:3.7.0`
- lftp - `docker build . --tag=lftp:1.0`
- cohortvcftozarr - `docker build . --tag=cohortvcftozarr:1.0`
- samplevcftozarr - `docker build . --tag=samplevcftozarr:1.0`
- sampleselectvariants - `docker build . --tag=sampleselectvarians:1.0`
- whatshap - `docker build . --tag=whatshap:1.0`
- shapeit4 - `docker build . --tag=shapeit4:4.2.1`
- bcftools - `docker build . --tag=bcftools:1.11`
- import - `docker build . --tag=sangerpathogens/irods:4.1.12`

