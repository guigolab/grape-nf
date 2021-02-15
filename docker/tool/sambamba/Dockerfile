# Dockerfile for sambamba, compiled from source
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER} AS builder

RUN apk update && \
    apk --no-cache add build-base \
                       curl \
                       git \
                       python3 \
                       libxml2 \
                       zlib-dev \
                       zlib-static

RUN apk --no-cache add -X http://dl-cdn.alpinelinux.org/alpine/edge/main \
                       -X http://dl-cdn.alpinelinux.org/alpine/edge/community \
                       ldc

ARG REVISION=v0.7.1

RUN git clone --recursive https://github.com/biod/sambamba.git --branch $REVISION \
    && cd sambamba && make static

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /sambamba/bin/sambamba-0.7.1 /usr/local/bin/sambamba
