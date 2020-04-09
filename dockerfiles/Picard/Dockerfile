FROM alpine:latest

ARG version=2.22.3

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

RUN apk update \
 && apk upgrade \
 && apk add openjdk8-jre-base \
 && apk add wget

RUN wget https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar
