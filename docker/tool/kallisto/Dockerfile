# Dockerfile for samtools
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER} AS builder

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"



RUN apk update && apk --no-cache add \
    gawk
