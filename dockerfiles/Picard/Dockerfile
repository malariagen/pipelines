FROM bashell/alpine-bash:latest

ARG version=2.22.3

WORKDIR /root

RUN apk update \
 && apk upgrade \
 && apk add openjdk8-jre-base \
 && apk add wget

RUN wget https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar
