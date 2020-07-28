# Mosquito SNP genotyping pipeline specification

* Version: 1.4.1
* Authors: Alistair Miles, Jim Stalker

This document specifies a pipeline for genotyping an individual sample
at a set of predefined SNP alleles, assuming the sample has already
been aligned against the appropriate reference sequence via the
[mosquito short read alignment
protocol](short-read-alignment-vector.md).


## Inputs

* Analysis-ready sequence read alignments for a single sample
  (produced by the [mosquito short read alignment
  pipeline](short-read-alignment-vector.md); BAM format)

* Alleles against which to genotype (VCF format)

* Reference sequence (FASTA format)


## Outputs

* A bgzipped VCF file with genotype calls for a single sample at the
  given alleles

* A zipped Zarr file converted from the VCF


## Software

* GATK version 3.7-0
* scikit-allel version 1.3.1
* Zarr version 2.4.0


## Pipeline


### Step: Genotyping

Genotype a single sample at given alleles using GATK UnifiedGenotyper
with the following parameters:

```bash
java -jar GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I {sample BAM} \
    --alleles {alleles VCF} \
    -R {reference sequence} \
    --out {output VCF} \
    --genotype_likelihoods_model BOTH \
    --genotyping_mode GENOTYPE_GIVEN_ALLELES \ 
    --heterozygosity 0.015 \
    --heterozygosity_stdev 0.05 \
    --indel_heterozygosity 0.001 \
    --downsampling_type BY_SAMPLE \
    -dcov 250 \
    --output_mode EMIT_ALL_SITES \
    --min_base_quality_score 17 \
    -stand_call_conf 0.0 \
    -contamination 0.0 \
    -A DepthPerAlleleBySample \
    -A RMSMappingQuality \
    -XA Coverage \
    -XA ExcessHet \
    -XA InbreedingCoeff \
    -XA MappingQualityZero \
    -XA HaplotypeScore \
    -XA SpanningDeletions \
    -XA FisherStrand \
    -XA StrandOddsRatio \
    -XA ChromosomeCounts \
    -XA BaseQualityRankSumTest \
    -XA MappingQualityRankSumTest \
    -XA QualByDepth \
    -XA ReadPosRankSumTest
```

Notes: 
 
1. Because we are genotyping a single sample against given alleles, we
   don't need many of the variant annotations that we would usually
   ask for when doing discovery. Reducing the number of annotations
   should reduce file size somewhat.

2. The '-contamination' option is set at zero (`0.0`) here rather than
   `0.05` used previously because of the fact that UnifiedGenotyper downsamples 
   reads carrying _non-ref_ alleles, leading to reference bias. The ref/alt depth 
   distribution for hets shifts from 50% to something lower. This can cause 
   systematic calling of genotypes as hom ref that would be called as het 
   without any downsampling.

3. When genotyping given alleles, GATK ignores any variant in the
   known variants VCF that is not marked as PASS in the FILTER
   field. However, we want to genotype at all variants in the known
   variants VCF, regardless of what is in the FILTER column. Therefore
   the known variants VCF must be created setting FILTER to PASS for
   all variants to get this to work as desired.

4. Given the sites VCF file is fixed for every sample, and we wish to generalise
   to future sets of sites/alleles, the VCF file describing sites and alleles 
   should be considered a parameter. This file for _A. gambiae_ (AgamP4) is available at 
   `gs://vo_agam_production/resources/observatory/ag.allsites.nonN.vcf.gz`
 
5. No further annotations or metrics need to be collected at this point.


### Step: VCF to Zarr conversion

Convert the VCF to Zarr format via the sample_vcf_to_zarr.py script,
with the following arguments:

* ``--sample {sample_identifier}``
* ``--field variants/MQ``
* ``--field calldata/GT``
* ``--field calldata/GQ``
* ``--field calldata/AD``

Additionally, for *Anopheles gambiae*, the ``--contig`` argument
should be provided once for each contig in the reference genome.

E.g.:

```bash
python sample_vcf_to_zarr.py \
    --input {path_to_vcf} \
    --output {path_to_zarr} \
    --sample {sample_identifier} \
    --field variants/MQ \
    --field calldata/GT \
    --field calldata/GQ \
    --field calldata/AD \
    --contig 2R \
    --contig 2L \
    --contig 3R \
    --contig 3L \
    --contig X \
    --contig Y_unplaced \
    --contig UNKN \
    --log {path_to_log} \
    --zip
```


## Implementation notes

* For all vector species we currently genotype at all sites in the
  genome and all possible SNP alleles, which means that the set of
  alleles provided to the genotyper comprises all genome sites where
  the reference allele is not "N" and all possible SNP alternate
  alleles are given in lexical order.


## Data structures

For the output of the "VCF to Zarr" conversion step, we expect that
this will create a Zarr hierarchy with the following structure:

* `/` [root group]
  * `{chromosome}` [group]
    * `variants` [group]
      * `MQ` [array; int8; shape (n_sites,)]
    * `calldata` [group]
      * `GT` [array; int8; shape (n_sites, 1, 2), dimensions correspond to (sites, samples, ploidy)]
      * `GQ` [array; int8; shape (n_sites, 1), dimensions correspond to (sites, samples)]
      * `AD` [array; int16; shape (n_sites, 1, 4), dimensions correspond to (sites, samples, alleles)]

There would be one `{chromosome}` group per chromosome, which for
*An. gambiae* would mean five groups named "2R", "2L", "3R", "3L" and
"X", corresponding to the five chromosomes in the AgamP4 reference
genome (excluding the "UNKN", "Y_unplaced" and mitochondrial genomes
which we don't genotype).


## Change log

* Version 1.4.1 - Modified the Zarr conversion step to use a
  script. Fixed missing MQ annotation in output VCF.

* Version 1.4.0 - Updated Zarr and scikit-allel versions. Corrected
  the fields parameters in the VCF to Zarr conversion (added
  variants/MQ).

* Version 1.3.0 - Version ported across from MalariaGEN vector-ops
  repo.
