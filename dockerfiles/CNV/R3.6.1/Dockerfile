## Platform tag ensures that images built on ARM-based machines (ex. M-series macs) are portable
FROM --platform="linux/amd64" ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Set parameteres
ARG R_VERSION="3.6.1"

WORKDIR /usr/local/

# Update and install dependencies and wet
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update \
    && apt-get install -y \
    build-essential \
    wget

# Install R and addtional dependencies
RUN apt-get build-dep -y r-base \
    && wget https://cloud.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz \
    && tar xzf R-${R_VERSION}.tar.gz \
    && cd "R-${R_VERSION}" \
    && ./configure --prefix=/opt/R/${R_VERSION}/ --enable-R-shlib --with-blas --with-lapack \
    && make \
    && make install

RUN rm -r R-${R_VERSION} \
    && rm R-${R_VERSION}.tar.gz

# Run an R script to install several libraries
COPY installRDeps.R /usr/local/
RUN /opt/R/${R_VERSION}/bin/Rscript /usr/local/installRDeps.R \
    && rm installRDeps.R

# Add R installation directory to path
ENV PATH="/opt/R/${R_VERSION}/bin/:$PATH"

# Copy source code
COPY . .



