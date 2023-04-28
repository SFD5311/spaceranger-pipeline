FROM continuumio/miniconda3

WORKDIR /opt


RUN apt-get update \
    && apt-get -y install curl \
    && apt-get -y install wget \
    && apt-get -y install git\
    && apt-get -y install \
    python3-pip \
 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

COPY data /opt

COPY docker/requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /out


COPY bin /opt
