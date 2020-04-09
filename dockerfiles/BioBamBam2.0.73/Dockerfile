FROM ubuntu:bionic-20191029
#Ubuntu builds are provided by repo. Could move to alpine but would need to build ourselves

RUN apt-get update \
 && apt-get install -y \
 tar \
 gzip \
 wget \
 xz-utils

WORKDIR /tmp

ENV file=biobambam2-2.0.73-release-20170620145717-x86_64-etch-linux-gnu.tar.gz

#Get biobambam2
RUN wget https://github.com/gt1/biobambam2/releases/download/2.0.73-release-20170620145717/${file}

RUN tar -C /tmp/ -xvzf /tmp/${file} \
&& rm ${file}

WORKDIR /tmp/biobambam2/2.0.73-release-20170620145717/x86_64-etch-linux-gnu/

RUN cp -r bin/* /usr/local/bin

RUN cp -r lib/* /usr/local/lib

WORKDIR /

RUN rm -r /tmp/biobambam2
