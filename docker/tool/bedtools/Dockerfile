# Dockerfile for BEDtools
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER} AS builder

ARG BEDTOOLS_VER=2.19.1

ENV _BEDTOOLS_VERSION ${BEDTOOLS_VER}

RUN apk update && apk --no-cache add \
    build-base \
    zlib-dev \
    curl \
    python \
    bash

RUN mkdir -p build/bedtools && \
    curl -L https://api.github.com/repos/arq5x/bedtools2/tarball/v$_BEDTOOLS_VERSION \
    | tar xzf - --strip-components=1 -C build/bedtools

RUN sed -i '111s/const/constexpr/' build/bedtools/src/utils/fileType/FileRecordTypeChecker.h && \
    cd build/bedtools && make

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /build/bedtools/bin/* /usr/local/bin/
