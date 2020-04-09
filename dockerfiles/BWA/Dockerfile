
FROM ubuntu:bionic-20191029

ARG version=v0.7.17

#Get GCC
RUN apt-get update \
    && apt-get install -y libc6 \
    zlib1g \
    build-essential \
    zlib1g-dev \
    make \
    git


WORKDIR /

RUN git clone https://github.com/lh3/bwa.git

WORKDIR /bwa

RUN git checkout tags/$version

RUN make

ENV PATH=$PATH:/bwa
