FROM bashell/alpine-bash@sha256:965a718a07c700a5204c77e391961edee37477634ce2f9cf652a8e4c2db858ff

ARG version=1.11

RUN apk update \
 && apk add --no-cache \
 build-base \
 zlib-dev \
 bzip2-dev \
 xz-dev \
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
  make install

# Install bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/${version}/bcftools-${version}.tar.bz2 && \
  tar xf bcftools-${version}.tar.bz2 && \
  rm bcftools-${version}.tar.bz2 && \
  cd bcftools-${version} && \
  ./configure --with-htslib=system && \
  make install
