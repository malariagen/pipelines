FROM bashell/alpine-bash@sha256:965a718a07c700a5204c77e391961edee37477634ce2f9cf652a8e4c2db858ff

ARG version=1.9

RUN apk update \
 && apk add --no-cache \
 build-base \
 zlib-dev \
 bzip2-dev \
 xz-dev \
 ncurses-dev \
 ca-certificates \
 wget

RUN wget -q https://github.com/samtools/samtools/releases/download/${version}/samtools-${version}.tar.bz2

RUN mkdir /samtools-${version}

RUN tar -C / -xjvf samtools-${version}.tar.bz2

WORKDIR /samtools-${version}/

RUN make

RUN mv /samtools-${version}/samtools /bin/ \
 && rm -rf /samtools*

RUN apk del build-base zlib-dev ca-certificates wget

RUN apk update \
 && apk upgrade \
 && apk add --no-cache python3
