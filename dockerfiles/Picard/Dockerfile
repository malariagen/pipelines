FROM bashell/alpine-bash@sha256:965a718a07c700a5204c77e391961edee37477634ce2f9cf652a8e4c2db858ff

ARG version=2.22.3

RUN apk update \
 && apk upgrade \
 && apk add openjdk8-jre-base \
 && apk add wget

WORKDIR /picard-${version}/

RUN wget https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar

WORKDIR /

RUN mv /picard-${version}/picard.jar /bin/ \
 && rm -rf /picard-${version}/
