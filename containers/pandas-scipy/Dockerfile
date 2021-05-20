FROM ubuntu:20.04

COPY requirements.txt /root/requirements.txt

RUN apt update && apt upgrade -y
RUN apt install -y python3 python3-pip

RUN pip3 install -r /root/requirements.txt
