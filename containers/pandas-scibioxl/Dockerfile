FROM ubuntu:20.04

ARG DEBIAN_FRONTEND="noninteractive"

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains pandas, biopython, openpyxl and scipy"
LABEL org.opencontainers.image.source https://github.com/catg-umag/ncov2019-ont-nf

COPY requirements.txt /root/requirements.txt

RUN apt-get update && apt upgrade -y \
	&& apt-get install -y python3 python3-pip unoconv \
	&& rm -fr /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir -r /root/requirements.txt
