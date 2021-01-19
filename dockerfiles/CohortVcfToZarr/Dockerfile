FROM python:3.7.2

RUN pip install --upgrade pip

COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY cohort_vcf_to_zarr.py .
COPY test_cohort_vcf_to_zarr.py .
COPY fixture ./fixture

