FROM ubuntu:bionic-20191029

RUN apt-get update \
  && apt-get install -y python3-pip \
  python3-dev \
  build-essential

RUN pip3 install --upgrade pip \
    && pip3 install whatshap \
    && export PATH=$HOME/.local/bin:$PATH
