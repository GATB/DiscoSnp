FROM debian:wheezy

MAINTAINER Pierre Peterlongo pierre.peterlongo@inria.fr

# Set discoSnp++ version
ENV MTG_VERSION 2.2.10

# Set noninteratve mode
ENV DEBIAN_FRONTEND noninteractive
ENV PACKAGES wget gcc g++ make cmake zlib1g-dev libboost-dev git

ENV DIR /opt
ENV SOURCEKS kissnp2
ENV BUILD build

WORKDIR ${DIR}

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends ${PACKAGES}

RUN git config --global http.sslVerify false

# clone the github repo
RUN git clone --recursive https://github.com/GATB/DiscoSnp.git

WORKDIR ${DIR}/${SOURCE}
RUN git submodule init

# Using an official release
RUN git checkout v${MTG_VERSION}
RUN git submodule update

RUN mkdir ${BUILD}
WORKDIR ${DIR}/${SOURCE}/${BUILD}

RUN cmake ..
RUN make

# symlink binary in /usr/local/bin
RUN ln -s ${DIR}/${SOURCE}/${BUILD}/bin/kissnp2 /usr/local/bin