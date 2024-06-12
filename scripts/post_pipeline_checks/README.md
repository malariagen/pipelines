# Post-pipeline checks 

## For Short-read alignment and SNP-genotyping pipeline output 

### check_missing_samples.sh
The script `check_missing_samples.sh` is used to quickly check whether all samples are present in an S3 bucket. The number of samples in the bucket and corresponding manifest are counted. Any missing sample is reported.

#### Running the script
To run the script, login to the farm as a user with permission to list files from S3 buckets via the `s3cmd` tool. Ensure that any manifest files are in the current working directory and named as `<batch>.tsv`. Then run:

```bash
./check_missing_samples.sh <batch_1> <batch_2> ...
```

where `<batch_N>` is the name of a production batch, e.g. 1177-VO-ML-LEHMANN-VMF00015.

### check_vcf_integrity.sh
The script `check_vcf_integrity.sh` is used to check the integrity of VCF files uploaded to S3. For a given batch, the script checks all uploaded vcf.gz files that are below a threshold filesize (2000000000 bytes = 2GB). It compares the number of lines corresponding to each chromosome in these VCF files to an expected number of lines. Any deviations are reported. 

#### Running the script
To run the script, login to the farm as a user with permission to list and download files from S3 buckets via the `s3cmd` tool. You can run the script with `bsub` to ensure that it has sufficient resources and survives as a long-running job (change the `-q` argument of `bsub` as required).

```bash
bsub -J <job_name> -o %J.o -e %J.e -M2000 -R'select[mem>2000] rusage[mem=2000]' -q long "./check_vcf_integrity.sh <batch>"
```

where `<batch>` is the name of the production batch, e.g. 1177-VO-ML-LEHMANN-VMF00015. Check the job output file `<job_id>.o` to view a summary of the number of lines for each chromosome and any deviations from expected values. To see only warnings/deviations, you can run `grep -E '(No rows|unexpected number)'` on script output.

