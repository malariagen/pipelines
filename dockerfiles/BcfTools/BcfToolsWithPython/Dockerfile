FROM python:3.7.2

ARG version=1.11

RUN apt-get update -qq -y && \
 apt-get install -y \
 build-essential \
 zlib1g-dev \
 bzip2 \
 xz-utils \
 ncurses-dev \
 ca-certificates \
 wget

WORKDIR /tmp

# Install htslib to get bgzip and tabix
RUN wget -q https://github.com/samtools/htslib/releases/download/${version}/htslib-${version}.tar.bz2 && \
  tar xf htslib-${version}.tar.bz2 && \
  rm htslib-${version}.tar.bz2 && \
  cd htslib-${version} && \
  ./configure && \
  make && \
  make install

# Install bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/${version}/bcftools-${version}.tar.bz2 && \
  tar xf bcftools-${version}.tar.bz2 && \
  rm bcftools-${version}.tar.bz2 && \
  cd bcftools-${version} && \
  ./configure --with-htslib=system && \
  make && \
  make install && \
  ldconfig
