FROM ubuntu:20.04

COPY VERSION VERSION

RUN apt update && apt upgrade -y
RUN apt install -y wget

RUN wget https://github.com/nextstrain/nextclade/releases/download/$(cat VERSION)/nextclade-Linux-x86_64 -O /usr/local/bin/nextclade
RUN chmod +x /usr/local/bin/nextclade

RUN mkdir /opt/nextclade_data
RUN cd /opt/nextclade_data && \
    wget -O - https://github.com/nextstrain/nextclade/archive/master.tar.gz | tar -xz --strip=2 "nextclade-master/data"
