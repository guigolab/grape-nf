# Dockerfile for RSeQC
#
FROM       grapenf/python

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ENV _RSEQC_VERSION 2.6.4

RUN apk add --no-cache --virtual=.build-dependencies build-base \
    zlib-dev \
    python-dev && \
    pip --no-cache-dir install rseqc==$_RSEQC_VERSION && \
    apk del .build-dependencies