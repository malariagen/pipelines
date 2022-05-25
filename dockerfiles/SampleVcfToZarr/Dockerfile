FROM python:3.8-slim-bullseye
ARG DEBIAN_FRONTEND=noninteractive

# Install OS dependencies
RUN apt-get update -qq -y && \
    apt-get install -y tabix

# Install Python dependencies
RUN pip3 install --upgrade pip
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Copy sample_vcf_to_zarr.py to WORKDIR
RUN mkdir /tools
WORKDIR /tools
COPY . .
