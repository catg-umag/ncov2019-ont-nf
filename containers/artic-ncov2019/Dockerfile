FROM mambaorg/micromamba:latest

ARG DEBIAN_FRONTEND="noninteractive"

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains ARTIC and nCoV-2019 primer schemes"
LABEL org.opencontainers.image.source https://github.com/catg-umag/ncov2019-ont-nf

SHELL ["/bin/bash", "-c"] 

USER root
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y wget procps && \
    rm -fr /var/lib/apt/lists/*

# download primer schemes for nCoV-2019
RUN mkdir /opt/artic-ncov2019
RUN cd /opt/artic-ncov2019 && \
    wget -qO- https://codeload.github.com/artic-network/artic-ncov2019/tar.gz/master \
        | tar -xz --strip=1 artic-ncov2019-master/primer_schemes/nCoV-2019

# install artic
RUN cat /etc/passwd
COPY env.yaml /tmp/env.yaml

RUN micromamba install -y -f /tmp/env.yaml && \
    micromamba clean -aqy

# reference files need to be indexed for the pipeline
USER root
RUN find /opt/artic-ncov2019/ -name '*.fasta' -exec /opt/conda/bin/samtools faidx {} \;

ENV PATH=/opt/conda/bin/:$PATH