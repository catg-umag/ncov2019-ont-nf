FROM python:3.11 as base

ARG DEBIAN_FRONTEND="noninteractive"

LABEL author="Diego Alvarez (dialvarezs@gmail.com)"
LABEL description="Contains pandas, biopython, openpyxl and scipy"
LABEL org.opencontainers.image.source https://github.com/catg-umag/ncov2019-ont-nf

COPY requirements.txt /tmp/requirements.txt

RUN apt-get update && apt upgrade

RUN pip3 install --no-cache-dir -r /tmp/requirements.txt


FROM python:3.11-slim as final

RUN apt-get update && apt upgrade -y \
	&& apt-get install -y \
		procps default-jre libreoffice-calc-nogui \
	&& rm -fr /var/lib/apt/lists/*

COPY --from=base /usr/local/ /usr/local
