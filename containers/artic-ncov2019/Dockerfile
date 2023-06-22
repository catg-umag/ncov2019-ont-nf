FROM mambaorg/micromamba:latest

ARG DEBIAN_FRONTEND="noninteractive"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains ARTIC and nCoV-2019 primer schemes"
LABEL org.opencontainers.image.source https://github.com/catg-umag/ncov2019-ont-nf

USER root

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y git procps && \
    rm -fr /var/lib/apt/lists/*

# download primer schemes
RUN cd /opt/ && \
    git clone https://github.com/artic-network/primer-schemes/ && \
    chown -R $MAMBA_USER:$MAMBA_USER primer-schemes

USER $MAMBA_USER

# install artic
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -yq -f /tmp/env.yaml && \
    micromamba clean -ayq

# reference files need to be indexed for the pipeline
RUN find /opt/primer-schemes/ -name '*.fasta' -exec samtools faidx {} \;

ENV PATH="/opt/conda/bin:${PATH}"
