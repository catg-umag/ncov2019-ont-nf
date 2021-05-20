FROM ubuntu:20.04

COPY requirements.txt /root/requirements.txt

RUN apt update && apt upgrade -y
RUN DEBIAN_FRONTEND="noninteractive" apt install -y python3 python3-pip unoconv

RUN pip3 install -r /root/requirements.txt
