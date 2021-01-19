# Mosquito phasing pipeline specification

* Version: 0.0.0
* Authors: Alistair Miles, Jon Brenas
* Status: DRAFT

This document specifies a pipeline for phasing SNP genotypes for a
cohort of multiple samples. The pipeline makes use of both read-backed
phasing and statistical (population) phasing to improve phasing
accuracy.


## Inputs

* A manifest specifying a set of samples to be phased.

* For each sample, a BAM file containing analysis-ready sequence read
  alignments, produced by the [mosquito short read alignment
  pipeline](short-read-alignment-vector.md).

* For each sample, a VCF file or zipped zarr file containing unphased
  genotypes, produced by the [mosquito SNP genotyping
  pipeline](snp-genotyping-vector.md). (@@TODO: decide which format
  will be more convenient to work from.)

* A file containing the set of sites and alleles to phase (VCF or Zarr
  format, @@TODO decide which). N.B., This will be a subset of the
  sites and alleles which were genotyped, because we only phase at
  sites that pass filtering, and we also can only phase biallelic
  variants due to limitations of the software.

* Reference genome sequence (FASTA format).

* Genetic map file.

* (OPTIONAL) A set of haplotypes previously phased for a different
  cohort of samples to use as a haplotype reference panel (VCF
  format).


## Outputs

* A bgzipped VCF file containing phased genotypes for all samples.

* A Zarr directory containing phased genotypes for all samples,
  converted from the VCF.


## Software

* [WhatsHap](https://whatshap.readthedocs.io/en/latest/) version 1.0

* [SHAPEIT4](https://odelaneau.github.io/shapeit4/) version 4.1.3

* TODO: other software needed to prepare and post-process files


## Pipeline

The pipeline comprises two main sub-pipelines, the first of which
performs read-backed phasing for each sample individually, and the
second of which performs simultaneous statistical phasing of all
samples together.


### Sub-pipeline: Read-backed phasing

This sub-pipeline can be run on each sample independently.


#### Step: Genotype data preparation

The input genotype data for each sample need to be subsetted down to
the set of variants at which phasing will be performed, which are
biallelic variants that pass filters. The input data may also contain
variants where only one alternate allele is observed, but where the
alternate allele observed is coded as ALT allele 2 or 3, and these
need to be recoded as biallelic genotypes (using 0 and 1 only).

This step takes as inputs:

* Genotypes for a single sample, as either VCF or zipped Zarr (@@TODO:
  decide which).

* A file specifying the subset of sites and alleles at which to phase,
  as Zarr.

This step then applies a subsetting and allele recoding operation to
the input genotypes, outputting a VCF file with biallelic genotypes
for a single sample via the sample_select_variants.py script.

E.g.:

```bash
python sample_select_variants.py \
      --sample-genotypes ~{sample_zarr} \
      --sites-called ~{called_sites_zarr} \
      --sites-selected ~{phased_sites_zarr} \
      --output ~{output_basename}.subset.vcf \
      --contig ~{contig} \
      --progress
```


#### Step: WhatsHap phase

Run WhatsHap ``phase`` on a single sample with the following command:

```bash
whatshap phase \
    -o sample_phased.vcf \
    --reference=reference.fasta \
    --internal-downsampling=@@TODO \
    sample_genotypes.vcf sample_alignments.bam
```

Notes:

* ``--internal-downsampling`` default is 15, docs say may be increased up to 20
  to improve phasing quality, but runtime increases exponentially, so
  some exploration of runtime for values in the range 15-20 might be
  useful to see what we can afford.

* A ``--chromosome`` option may be given to phase each chromosome
  separately, which may be useful for parallelisation.


#### Step: WhatsHap stats

Run WhatsHap ``stats`` to compute phasing statistics for a phased
sample:

```bash
whatshap stats \
    --chr-lengths=CHR_LENGTHS \
    --tsv=stats.tsv \
    --gtf=blocks.gtf \
    sample_phased.vcf
```

Notes:

* A ``--chromosome`` option may be given to run each chromosome
  separately, which may be useful for parallelisation.


### Sub-pipeline: Statistical phasing

This sub-pipeline merges the partially-phased genotypes for each
sample produced by WhatsHap into a single VCF and performs statistical
(population) phasing using SHAPEIT4.


#### Step: Merge VCFs

Merge the outputs from the WhatsHap phase step into a single
multi-sample VCF.

@@TODO which tool/command to use.


#### Step: SHAPEIT4

Run SHAPEIT4 on the merged VCF containing partially phased genotypes
from WhatsHap.

```bash
shapeit4 \
    --input merged.vcf \
    --map GMAP \
    --region REGION \
    --window WINDOW \
    --thread N_THREADS \
    --mcmc-iterations MCMC_ITERATIONS \
    --pbwt-depth PBWT_DEPTH \
    --window WINDOW \
    --sequencing \
    --reference REFERENCE.vcf \
    --use-PS 0.0001 \
    --log phased.log \
    --output phased.vcf.gz
```

Notes:

* SHAPEIT4 can be run on a specific genome region (``--region``),
  allowing some parallelisation by genome region. @@TODO discuss
  appropriate region size.

* A genetic map file is needed (``--map``). We have one for
  *An. gambiae*. For other species where no map is available a map can
  *be created assuming constant recombination rate.

* SHAPEIT4 can use multi-threading (``--threads``) which may be
  important especially when phasing larger cohorts.

* SHAPEIT4 allows tuning of the MCMC iteration scheme
  (``--mcmc-iterations``), which allows the tradeoff between accuracy
  and speed to be explored, see
  [documentation](https://odelaneau.github.io/shapeit4/#documentation)
  point 5. The default is ``5b,1p,1b,1p,1b,1p,5m`` but an alternative
  offering greater accuracy is suggested
  ``10b,1p,1b,1p,1b,1p,1b,1p,10m``. @@TODO spike to get a sense of
  performance cost of increasing number of burn-in and/or main
  iterations.

* SHAPEIT4 allows tuning of the number of conditioning haplotypes via
  the ``--pbwt-depth`` parameter, with default 4 and higher values
  providing more accuracy at slower performance. @@TODO spike
  performance cost of increasing.

* SHAPEIT4 allows tuning of the size of the window used for phasing
  (``--window``) which is given in cM, default 2.5. @@TODO check this
  is sensible for *An. gambiae*.

* Parameters may need to be tuned to adjust for high variant
  density. The option ``--sequencing`` alters default parameters for
  PBWT and IBD2 mapping to be appropriate for sequencing data, and so
  should be a good baseline to start from, but some tuning of
  individual parameters may still be needed. @@TODO figure out how
  --sequencing changes defaults for other params.

* The ``--reference`` parameter is optional and given only if phasing
  against a reference panel.

* The parameter ``--use-PS`` tells SHAPEIT4 to use the phasing
  information provided by WhatsHap.


#### Step: Ligate regions?

@@TODO If SHAPEIT4 is run in separate genome regions, do we need to
ligate the regions back together? This was previously needed for
SHAPEIT2 but I can't see anything about this on the SHAPEIT4 docs. Or
should we instead phase entire chromosomes at a time?

Output should be a whole-genome VCF, or one VCF per chromosome also
acceptable.


#### Step: VCF to Zarr

For ease of downstream analysis, convert the final VCF to Zarr via the
scikit-allel vcf_to_zarr function.

```
@@TODO function params
```

N.B., there is no need to zip the outputs of this step, we will use
them unzipped.


## Implementation notes

@@TODO something about chunking and parallelisation

* SHAPEIT4 has optimisations to use AVX2 instruction set, and so
  compiling and running the program on CPUs supporting AVX2 is
  preferable. See also [SHAPEIT4 installation
  notes](https://odelaneau.github.io/shapeit4/#installation).


## Changes

* Version 0.0.0 - First draft.
