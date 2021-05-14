FROM ubuntu:20.04

RUN apt update && apt upgrade -y
RUN apt install -y wget

RUN mkdir /opt/artic-ncov2019
RUN cd /opt/artic-ncov2019 && \
    wget -qO- https://codeload.github.com/artic-network/artic-ncov2019/tar.gz/master \
        | tar -xz --strip=1 artic-ncov2019-master

ENV MAMBA_ROOT_PREFIX=/opt/micromamba/
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest \
        | tar -xvj bin/micromamba
RUN micromamba env create -f /opt/artic-ncov2019/environment.yml
RUN micromamba clean -aqy

ENV PATH=/opt/micromamba/envs/artic-ncov2019/bin:$PATH