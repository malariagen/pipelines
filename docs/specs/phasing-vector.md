# Phasing (vector) pipeline specification

* Version: 0.0.0
* Authors: Alistair Miles, Jon Brenas
* Status: WORK IN PROGRESS

This document specifies a pipeline for phasing SNP genotypes for a set
of samples. The pipeline makes use of both read-backed phasing and
statistical (population) phasing to improve phasing accuracy.


## Inputs

* A manifest specifying a set of samples to be phased.

* For each sample, a BAM file containing analysis-ready sequence read
  alignments, produced by the short read alignment pipeline.

* For each sample, a VCF file or zipped zarr file containing unphased
  genotypes, produced by the SNP genotyping pipeline. @@TODO: decide
  which format will be more convenient to work from.

* A file containing the set of sites and alleles to phase (VCF or Zarr
  format, @@TODO decide which). (This will be a subset of the sites
  and alleles which were genotyped, because we only phase at sites
  that pass filtering, and we also can only phase biallelic variants
  due to limitations of the software.)

* Reference genome sequence (FASTA format).

* Genetic map file.

* (OPTIONAL) A set of haplotypes previously phased for a different set
  of samples to use as a haplotype reference panel (VCF format).


## Outputs

* A bgzipped VCF file containing phased genotypes for all samples.

* A Zarr directory containing phased genotypes, converted from the VCF.


## Software

* [WhatsHap](https://whatshap.readthedocs.io/en/latest/) version 0.18

* [SHAPEIT4](https://odelaneau.github.io/shapeit4/) version 4.1.3

* TODO: other software needed to prepare and post-process files


## Pipeline

The pipeline comprises two main sub-pipelines, the first of which
performs read-backed phasing for each sample individually, and the
second of which performs simultaneous phasing of all samples together.


### Sub-pipeline: Read-backed phasing

This sub-pipeline can be run on each sample independently.


#### Step: Genotype data preparation

The input genotype data for each sample need to be subsetted down to
the set of variants at which phasing will be performed. The input data
will also contain multiallelic variants, and these need to be recoded
as biallelic genotypes (using 0 and 1 only).

This step takes as inputs:

* Genotypes for a single sample, as either VCF or zipped Zarr (@@TODO:
  decide which).

* A file specifying the subset of sites and alleles at which to phase,
  as either VCF or Zarr (@@TODO: decide which).

This step then applies a subsetting and allele recoding operation to
the input genotypes, outputting a VCF file with biallelic genotypes
for a single sample.

@@TODO figure out which tool and command to use


#### Step: WhatsHap

@@TODO


### Sub-pipeline: Statistical phasing

This sub-pipeline merges the partially-phased genotypes for each
sample into a single VCF and performs statistical (population) phasing
using SHAPEIT4.


#### Step: Merge VCFs

@@TODO


#### Step: SHAPEIT4

@@TODO


#### Step: VCF to Zarr

@@TODO


## Implementation notes

@@TODO something about chunking and parallelisation


## Changes

* 