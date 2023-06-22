FROM mambaorg/micromamba:latest

ARG DEBIAN_FRONTEND="noninteractive"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains ARTIC and nCoV-2019 primer schemes"
LABEL org.opencontainers.image.source https://github.com/catg-umag/ncov2019-ont-nf

USER root

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y wget procps && \
    rm -fr /var/lib/apt/lists/*

RUN mkdir /opt/artic-ncov2019 && \
    chown $MAMBA_USER:$MAMBA_USER /opt/artic-ncov2019

USER $MAMBA_USER

# download primer schemes for nCoV-2019
RUN cd /opt/artic-ncov2019 && \
    wget -qO- https://codeload.github.com/artic-network/artic-ncov2019/tar.gz/master \
        | tar -xz --strip=1 artic-ncov2019-master/primer_schemes/nCoV-2019

# install artic
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

RUN micromamba install -yq -f /tmp/env.yaml && \
    micromamba clean -aqy

# reference files need to be indexed for the pipeline
RUN find /opt/artic-ncov2019/ -name '*.fasta' -exec samtools faidx {} \;

ENV PATH="/opt/conda/bin:${PATH}"
