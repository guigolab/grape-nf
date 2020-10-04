# Dockerfile for bamstats
#
ARG ALPINE_VER=3.11
FROM golang:alpine${ALPINE_VER} AS builder

ARG BAMSTATS_VER=0.3.4

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ENV _BAMSTATS_VERSION ${BAMSTATS_VER}

RUN apk update && \
    apk --no-cache add git

RUN GO111MODULE=on go get \
    -ldflags "-X github.com/guigolab/bamstats.PreVersionString=" \
    github.com/guigolab/bamstats/cmd/bamstats@v${_BAMSTATS_VERSION}

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /go/bin/bamstats /usr/local/bin/
