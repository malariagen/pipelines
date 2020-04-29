
# Dockerfile images

### BWA
Supports the `build-arg` `version` that represents the tag name found in the repo
https://github.com/lh3/bwa. Default version v0.7.17.

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
- For GATK 4.X use https://hub.docker.com/r/broadinstitute/gatk/

## Mosquito short read alignment pipeline

- bwa 0.7.15 - `docker build . --tag=bwa:0.7.15 --build-arg version=v0.7.15`
- samtools 1.4.1 (htslib 1.4.1) - `docker build . --tag=samtools:1.4.1 --build-arg version=1.4.1`
- picard 2.9.2 - ` docker build . --tag=picard:2.9.2 --build-arg version=2.9.2`
- biobambam 2.0.73 - `docker build . --tag=biobambam:2.0.73`
- GATK 3.7-0 - `docker pull broadinstitute/gatk3:3.7-0`
