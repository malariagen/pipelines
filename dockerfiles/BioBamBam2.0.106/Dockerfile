FROM ubuntu:bionic-20191029
#Ubuntu builds are provided by repo. Could move to alpine but would need to build ourselves

RUN apt-get update \
 && apt-get install -y \
 tar \
 gzip \
 wget \
 xz-utils

WORKDIR /tmp

ENV file=biobambam2_x86_64-linux-gnu_2.0.106.tar.xz

#Get biobambam2
RUN wget https://gitlab.com/german.tischler/biobambam2/uploads/8006c668e1fad2875f1fd2d952a502fb/${file}

RUN tar -C /tmp/ -vxJf /tmp/${file} \
&& rm ${file}

WORKDIR /tmp/biobambam2/x86_64-linux-gnu/2.0.106

RUN cp -r bin/* /usr/local/bin

RUN cp -r lib/* /usr/local/lib

WORKDIR /

RUN rm -r /tmp/biobambam2
