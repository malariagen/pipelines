# SNP genotyping (vector) pipeline specification

* Version: 1.4.0
* Authors: Alistair Miles, Jim Stalker

This document specifies a pipeline for genotyping an individual sample
at a set of predefined SNP alleles, assuming the sample has already
been aligned against the appropriate reference sequence via the short
read alignment protocol.


## Inputs

* Analysis-ready sequence read alignments for a single sample
  (produced by the short read alignment pipeline; BAM format)

* Alleles against which to genotype (VCF format)

* Reference sequence (FASTA format)


## Outputs

* A bgzipped VCF file with genotype calls for a single sample at the
  given alleles

* A zipped Zarr file converted from the VCF


## Software

* GATK version 3.7-0
* scikit-allel version 1.2.1
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
    -XA RMSMappingQuality \
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

2. I've set the '-contamination' option to zero (0.0) here rather than
   0.05 we've used previously because I didn't realise that this
   option down-samples only alternate alleles (I thought it
   downsampled both alt and ref alleles), and so setting contamination
   > 0 will create bias towards reference allele. We will deal with
   contamination separately.

3. When genotyping given alleles, GATK ignores any variant in the
   known variants VCF that is not marked as PASS in the FILTER
   field. However, we want to genotype at all variants in the known
   variants VCF, regardless of what is in the FILTER column. Therefore
   the known variants VCF must be created setting FILTER to PASS for
   all variants to get this to work as desired.


### Step: VCF to Zarr conversion

Convert the VCF to Zarr format via the scikit-allel
[vcf_to_zarr()](https://scikit-allel.readthedocs.io/en/stable/io.html#allel.vcf_to_zarr)
function. Note the following implementation details:

* Only the ``variants/MQ``, ``calldata/GT``, ``calldata/GQ`` and
  ``calldata/AD`` fields are required in the output.

* Data types should be as follows: ``variants/MQ``: i1 (single byte
  signed integer), ``calldata/GT``: i1; ``calldata/GQ``: i1;
  ``calldata/AD``: i2

* Compressor should be ``zarr.Blosc(cname='zstd', clevel=5,
  shuffle=-1)``

* Chunk length should be ``2**18``

* ``alt_number`` parameter should be 3 (fix number of alt alleles)

* Data should be grouped by contig (chromosome). This means calling
  ``vcf_to_zarr()`` once for each contig in the reference genome,
  specifying the name of the contig as the value of the ``region`` and
  ``group`` parameters.

Here is an example of the appropriate function call in Python for a
single contig (this should be called for all contigs in the genome):

```python
import sys
import allel
import zarr


def build_sample_zarr(input_path, output_path, contig):
    allel.vcf_to_zarr(
        input=input_path,
        output=output_path,
        group=contig,
        region=contig,
        compressor=zarr.Blosc(cname='zstd', clevel=5, shuffle=0),
        overwrite=True,
        fields=['calldata/GT', 'calldata/GQ', 'calldata/AD', 'variants/MQ'],
        types={'calldata/GT': 'i1',
               'calldata/GQ': 'i1',
               'calldata/AD': 'i2',
	       'variants/MQ': 'i1'},
        alt_number=3,
        chunk_length=2**18,
        log=sys.stdout,
    )

```


### Step: Zip Zarr

Once the zarr conversion is complete, create a zip archive from the
contents of the root folder of the zarr directory. E.g., if
``{sample}`` is the sample ID and ``{sample}.zarr`` is the name of the
output directory used in calls to ``vcf_to_zarr()``, then do:

```bash
cd {sample}.zarr
zip -rmT0 {sample}.zarr.zip .
```


## Implementation notes

* For all vector species we currently genotype at all sites in the
  genome and all possible SNP alleles, which means that the set of
  alleles provided to the genotyper comprises all genome sites where
  the reference allele is not "N" and all possible SNP alternate
  alleles are given in lexical order.


## Appendix: Data structures

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

* Version 1.4.0 - Updated Zarr and scikit-allel versions. Corrected
  the fields parameters in the VCF to Zarr conversion (added
  variants/MQ).

* Version 1.3.0 - Version ported across from MalariaGEN vector-ops
  repo.
