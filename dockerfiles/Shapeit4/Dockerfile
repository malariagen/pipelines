FROM ubuntu:20.04
# Suppress interactive output from installs

ENV DEBIAN_FRONTEND=noninteractive
ARG shapeit4_version=4.2.1
ARG htslib_version=1.11
ARG boost_version_with_dots=1.74.0
ARG boost_version=1_74_0

RUN apt-get update && \
    apt-get -y --no-install-recommends install \
      ca-certificates \
      libbz2-dev \
      libcurl4-openssl-dev \
      liblzma-dev \
      libssl-dev \
      libncurses5-dev \
      autoconf \
      automake \
      bzip2 \
      gcc \
      g++ \
      make \
      wget \
      zlib1g-dev \
      libboost-all-dev

WORKDIR /tmp
# Install htslib to get bgzip and tabix
RUN wget -q https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 && \
  tar xf htslib-${htslib_version}.tar.bz2 && \
  rm htslib-${htslib_version}.tar.bz2 && \
  cd htslib-${htslib_version} && \
  ./configure && \
  make install

RUN wget -q https://boostorg.jfrog.io/artifactory/main/release/${boost_version_with_dots}/source/boost_${boost_version}.tar.gz && \
  tar xf boost_${boost_version}.tar.gz && \
  rm boost_${boost_version}.tar.gz && \
  cd boost_${boost_version} && \
  ./bootstrap.sh --with-libraries=iostreams,program_options && \
  ./b2 install

RUN wget -q https://github.com/odelaneau/shapeit4/archive/v${shapeit4_version}.tar.gz && \
  tar xf v${shapeit4_version}.tar.gz && \
  rm v${shapeit4_version}.tar.gz && \
  cd shapeit4-${shapeit4_version} && \
  sed -i s^HTSLIB_INC=\$\(HOME\)/Tools/htslib-1.11^HTSLIB_INC=/usr/local/include/htslib^ makefile && \
  sed -i s^HTSLIB_LIB=\$\(HOME\)/Tools/htslib-1.11/libhts.a^HTSLIB_LIB=/usr/local/lib/libhts.a^ makefile && \
  sed -i s^BOOST_INC=/usr/include^BOOST_INC=/usr/local/include/boost/^ makefile && \
  sed -i s^BOOST_INC=/usr/include^BOOST_INC=/usr/local/include/boost/^ makefile && \
  sed -i s^BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a^BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a^ makefile && \
  sed -i s^BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a^BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a^ makefile && \
  make && \
  cp /tmp/shapeit4-${shapeit4_version}/bin/shapeit4.2 /usr/local/bin/shapeit4

WORKDIR /usr/local/bin

