# Mosquito CNV pipeline specification

* Version: 1.0.0
* Authors: Alistair Miles, Eric Lucas, Kevin Palis, Jessica Way

This document specifies a pipeline for calling Copy Number Variants (CNVs) for a
cohort of multiple samples. The pipeline inclludes an HMM step followed by coverage calls and target regions pipelines to improve acuracy. 


## Inputs

* A manifest specifying a set of samples to be called.

* For each sample, a BAM file containing analysis-ready sequence read
  alignments, produced by the [mosquito short read alignment
  pipeline](short-read-alignment-vector.md).

* A  unique sample set ID

* Reference genome sequence (FASTA format).


## Outputs

* CNV HMM output

* CNV Coverage Calls output

* CNV targeted regions output


## Software

* GATK version 3.7-0


## Pipeline

### Sub-pipeline: CNV HMM 

#### Step: Windowed Coverage
Software: Python virtualenv cnv37
Steps in CNV_pipeline/scripts/get_windowed_coverage_and_diagnostic_reads.sh
```bash

```

#### Step: Calculate Coverage Summary Stats
Software: Python virtualenv cnv37
Steps in CNV_pipeline/scripts/get_coverage_stats_by_sample_set_vobs.sh
```bash

```

#### Step: CNV HMM
Software: Python virtualenv cnv37
Steps in CNV_pipeline/scripts/coverage_HMM_vobs.sh
```bash

```

#### Step: Join Coverage Variance Files
Steps in CNV_scripts/scripts/join_species_coverage_variance_files_vobs.sh
```bash

```


### Sub-pipeline: CNV Coverage Calls

#### Step: CNV Coverage Calls
Software: R version 3.6.1
Steps in CNV_pipeline/scripts/coverage_CNVs_vobs.sh
```bash

```


### Sub-pipeline: CNV Targeted Regions 

#### Step: Extract Diagnostic Reads
Note: this is currently combined with the windowed coverage step at the beginning of the pipeline

Software: Python virtualenv cnv37
Steps in CNV_pipeline/scripts/get_windowed_coverage_and_diagnostic_reads.sh
```bash

```

#### Step: Target Regions CNV Calls
Software: R version 3.6.1
Steps in CNV_pipeline/scripts/target_regions_analysis_vobs.sh
```bash

```

#### Step: Join CNV Output Files
Steps in CNV_scripts/scripts/join_CNV_output_files_vobs.sh
```bash

```

#### Step: Run Modal CNVs
Software: R version not specified
Steps in CNV_pipeline/scripts/modal_CNVs.sh
```bash

```
**Notes:**
* Called separately for each species
* Not using the "vobs" version of the script, is this correct?

## Implementation notes



## Change log

