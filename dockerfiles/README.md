
# Dockerfile images

### BWA
Supports the `build-arg` `version` that represents the tag name found in the BWA repo
https://github.com/lh3/bwa. Defaults to v0.7.17.

**Currently only builds for v0.7.16 and above**

### Samtools
Supports the `build-arg` `version` that represents the tag name found in the BWA repo
https://github.com/samtools/samtools. Defaults to 1.9.

### BioBamBam

### Picard

### GATK

## Mosquito short read alignment pipeline

- bwa 0.7.15 - unsupported
- samtools 1.4.1 (htslib 1.4.1) - `docker build . --tag=samtools1.4.1 --build-arg version=1.4.1`
- picard 2.9.2 - broad tag history dose not go back this far!
- biobambam 2.0.73
- GATK 3.7-0
