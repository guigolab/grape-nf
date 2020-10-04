# Dockerfile for the grape-nf base image
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER}

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache \
    bash \
    bzip2 \
    ca-certificates \
    coreutils \
    curl \
    findutils \
    libgomp \
    libstdc++ \
    pigz \
    procps \
    tar \
    xz-libs
