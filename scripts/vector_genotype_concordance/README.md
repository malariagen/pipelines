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


## License
vector_genotype_concordance is free software, licensed under [MIT](https://github.com/malariagen/pipelines/blob/master/LICENSE).

