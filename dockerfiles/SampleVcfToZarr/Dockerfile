FROM python:3.7.2
ARG DEBIAN_FRONTEND=noninteractive

ARG HTSLIB_VERSION=1.11
ARG HTSLIB_URL="https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz"
ARG HTSLIB_DIR="/opt/htslib"
ARG BUILD_TOOLS="build-essential autoconf"

# Install htslib dependencies
RUN apt-get update -qq -y && \
    apt-get purge -y libssl-dev && \
    apt-get install -y ${BUILD_TOOLS} \
                       libssl1.0-dev \
                       curl libcurl4-openssl-dev

# Install htslib
RUN mkdir -p ${HTSLIB_DIR} && \
    cd ${HTSLIB_DIR} && \
    curl -fsSL ${HTSLIB_URL} | tar xzf - --strip-components=1 && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install

# Clean up
RUN apt-get remove -y ${BUILD_TOOLS} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install sample_vcf_to_zarr.py dependencies
# Note that numpy is installed separately to avoid a wheel build failure for scikit-allel
RUN pip install --upgrade pip \
 && pip3 install numpy \
 && pip3 install \
      scikit-allel==1.3.1 \
      zarr==2.4.0

# Copy sample_vcf_to_zarr.py to WORKDIR
RUN mkdir /tools
WORKDIR /tools
COPY . .
