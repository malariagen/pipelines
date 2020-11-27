FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
ARG IRODS_VERSION=4.2.8
RUN apt-get update \
    && apt-get install -yq \
      software-properties-common \
      wget \
      python3 \
      python3-pip \
    && pip3 install pandas \
    && wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
    && echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" | tee /etc/apt/sources.list.d/renci-irods.list \
    && apt-get update \
    && apt-get install -yq irods-icommands=${IRODS_VERSION} \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 10 \
