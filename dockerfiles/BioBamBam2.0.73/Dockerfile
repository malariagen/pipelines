FROM ubuntu:bionic-20191029

RUN apt-get update \
 && apt-get install -y \
 tar \
 gzip \
 wget \
 xz-utils

WORKDIR /tmp

#Get biobambam2
RUN wget https://github.com/gt1/biobambam2/releases/download/2.0.73-release-20170620145717/biobambam2-2.0.73-release-20170620145717-x86_64-etch-linux-gnu.tar.gz

RUN ls

RUN tar -C /tmp/ -xvzf /tmp/biobambam2-2.0.73-release-20170620145717-x86_64-etch-linux-gnu.tar.gz \
&& rm biobambam2-2.0.73-release-20170620145717-x86_64-etch-linux-gnu.tar.gz

WORKDIR /tmp/biobambam2/2.0.73-release-20170620145717/x86_64-etch-linux-gnu/

RUN ls

RUN cp -r bin/* /usr/local/bin

RUN cp -r lib/* /usr/local/lib

WORKDIR /

RUN rm -r /tmp/biobambam2
