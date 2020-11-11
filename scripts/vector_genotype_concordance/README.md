# Vector genotype concordance
A tool to compute and summarize genotype concordance for malaria vectors

## Pre-requiste
This package depends on scikit-allel which in turns requires the python development libraries to be installed. On ubuntu, run
```
sudo apt-get install python3-dev
```

## Install
Clone the repository locally and ideally install in a virtual environment:
```
# Clone the repo locally
git clone https://github.com/malariagen/pipelines
cd scripts/vector_genotype_concordance

# Create a virtual environment, if using virtualenv
virtualenv venv

# Create a virtual enviromnet, using venv:
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install vector_genotype_concordance
pip install wheel numpy
pip install .

# Run the tests
python setup.py nosetests
```

## Uninstall
```
pip uninstall vector_genotype_concordance
```

## Usage
The application has two modes: count and summarize.  Count will categorize genotype concordance, summarize produces concordance summary calculation at sample, category, chromosome and global level.
To get help:
```
vector_genotype_concordance.py -h
```
Count usage:
```
vector_genotype_concordance.py count -c '/path/to/zarrs/{sample}/{sample}.genotypes.zarr.zip' \
                                     -t '/path/to/zarrs/{sample}.zarr.zip'  \
                                     -s AA0052-C AB0252-C AC0010-C AJ0037-C AN0131-C AN0280-Cx AN0326-C AR0078-C \
                                     -o /some/dir/results.tsv 
```
Summary usage:
```
vector_genotype_concordance.py summarize -i results.*.tsv -o /some/dir/10_samples
```

## License
vector_genotype_concordance is free software, licensed under [MIT](https://github.com/malariagen/pipelines/blob/master/LICENSE).

