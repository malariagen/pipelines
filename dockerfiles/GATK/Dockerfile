FROM broadinstitute/gatk3:3.7-0
ARG DEBIAN_FRONTEND=noninteractive

ARG HTSLIB_VERSION=1.11
ARG HTSLIB_URL="https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz"
ARG HTSLIB_DIR="/opt/htslib"
ARG BUILD_TOOLS="build-essential autoconf"

# Change repositories that apt uses (to cope with Debian Jessie in base image)
RUN echo "deb http://cdn-fastly.deb.debian.org/debian/ jessie main" > /etc/apt/sources.list && \
    echo "deb-src http://cdn-fastly.deb.debian.org/debian/ jessie main" >> /etc/apt/sources.list && \
    echo "deb http://security.debian.org/ jessie/updates main" >> /etc/apt/sources.list && \
    echo "deb-src http://security.debian.org/ jessie/updates main" >> /etc/apt/sources.list && \
    echo "deb [check-valid-until=no] http://archive.debian.org/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list && \
    echo "deb-src [check-valid-until=no] http://archive.debian.org/debian jessie-backports main" >> /etc/apt/sources.list.d/jessie-backports.list 

# Install dependencies
RUN apt-get update -o Acquire::Check-Valid-Until=false -qq -y && \
    apt-get install -y ${BUILD_TOOLS} \
                       libcurl4-openssl-dev libssl-dev

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
