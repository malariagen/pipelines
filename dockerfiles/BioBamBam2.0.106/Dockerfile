FROM ubuntu:bionic-20191029

RUN apt-get update \
 && apt-get install -y \
 tar \
 gzip \
 wget \
 xz-utils

WORKDIR /tmp

#Get biobambam2
RUN wget https://gitlab.com/german.tischler/biobambam2/uploads/8006c668e1fad2875f1fd2d952a502fb/biobambam2_x86_64-linux-gnu_2.0.106.tar.xz

RUN tar -C /tmp/ -vxJf /tmp/biobambam2_x86_64-linux-gnu_2.0.106.tar.xz \
&& rm biobambam2_x86_64-linux-gnu_2.0.106.tar.xz

WORKDIR /tmp/biobambam2/x86_64-linux-gnu/2.0.106

RUN cp -r bin/* /usr/local/bin

RUN cp -r lib/* /usr/local/lib

WORKDIR /

RUN rm -r /tmp/biobambam2
