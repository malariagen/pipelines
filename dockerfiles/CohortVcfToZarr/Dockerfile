FROM python:3.7.2

RUN pip install --upgrade pip

COPY requirements.txt .
RUN pip3 install numpy==1.19.1
RUN pip3 install zarr==2.4.0
RUN pip3 install scikit-allel==1.3.1

RUN mkdir /tools
WORKDIR /tools

COPY cohort_vcf_to_zarr.py .
COPY test_cohort_vcf_to_zarr.py .
COPY fixture ./fixture

