FROM ubuntu:20.04

COPY VERSION /root/VERSION

RUN apt update && apt upgrade -y
RUN apt install -y wget

# download primer schemes for nCoV-2019
RUN mkdir /opt/artic-ncov2019
RUN cd /opt/artic-ncov2019 && \
    wget -qO- https://codeload.github.com/artic-network/artic-ncov2019/tar.gz/master \
        | tar -xz --strip=1 artic-ncov2019-master/primer_schemes/nCoV-2019

# setup micromamba
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest \
        | tar -xvj bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba/
RUN mkdir $MAMBA_ROOT_PREFIX
ENV PATH=/opt/micromamba/bin/:$PATH

# install artic
RUN eval "$(micromamba shell hook -s posix)" && \
    micromamba activate && \
    micromamba install artic=$(cat /root/VERSION) -c bioconda -c conda-forge
RUN micromamba clean -aqy

# reference files need to be indexed for the pipeline
RUN find /opt/artic-ncov2019/ -name '*.fasta' -exec samtools faidx {} \;