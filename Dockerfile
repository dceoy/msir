FROM dceoy/samtools:latest AS builder

ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD . /tmp/msir

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        python3.7 python3.7-distutils

RUN set -e \
      && /usr/bin/python3.7 /tmp/get-pip.py \
      && pip install -U --no-cache-dir pip /tmp/msir


FROM ubuntu:latest

ENV DEBIAN_FRONTEND noninteractive
COPY --from=builder /usr/local /usr/local

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        python3.7 \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["msir"]
