# Dockerfile for RSeQC
#
FROM grapenf/python

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ARG RSEQC_VER=2.6.4

ENV _RSEQC_VERSION ${RSEQC_VER}

RUN apk add --no-cache --virtual=.build-dependencies build-base \
    zlib-dev \
    python-dev && \
    pip --no-cache-dir install rseqc==$_RSEQC_VERSION && \
    apk del .build-dependencies