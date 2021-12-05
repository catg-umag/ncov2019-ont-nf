FROM ubuntu:20.04

ARG DEBIAN_FRONTEND="noninteractive"

COPY VERSION /root/VERSION

RUN apt-get update && apt upgrade -y \
    && apt-get install -y wget \
    && rm -fr /var/lib/apt/lists/*

RUN wget https://github.com/nextstrain/nextclade/releases/download/$(cat /root/VERSION)/nextclade-Linux-x86_64 \
        -O /usr/local/bin/nextclade
RUN chmod +x /usr/local/bin/nextclade