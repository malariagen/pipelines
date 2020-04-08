#FROM bashell/alpine-bash:latest
FROM alpine:latest

ARG version=v0.7.17

# make sure the package repository is up to date
RUN apk update \
 && apk upgrade \
 && apk add bash \
 && rm -rf /var/cache/*/* \
 && echo "" > /root/.ash_history

# change default shell from ash to bash
RUN sed -i -e "s/bin\/ash/bin\/bash/" /etc/passwd

ENV LC_ALL=en_US.UTF-8

WORKDIR /root

#Get GCC & git
RUN apk update \
 && apk upgrade \
 && apk add build-base \
 && apk add gcc \
 && apk add git

#Get required packages
RUN apk update \
 && apk upgrade \
 && apk add libxml2-dev \
 && apk add gd-dev \
 && apk add gzip

WORKDIR /

RUN git clone https://github.com/lh3/bwa.git

WORKDIR /bwa

RUN git checkout tags/$version

RUN make

ENV PATH=$PATH:/bwa
