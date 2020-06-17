# Dockerfile for samtools
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER} AS builder

ARG SAMTOOLS_VER=1.3.1

ENV _SAMTOOLS_VERSION ${SAMTOOLS_VER}

RUN apk update && apk --no-cache add \
    build-base \
    zlib-dev \
    ncurses-dev \
    curl

RUN mkdir -p build/samtools && \
    curl -L https://github.com/samtools/samtools/releases/download/$_SAMTOOLS_VERSION/samtools-$_SAMTOOLS_VERSION.tar.bz2 \
    | tar xj --strip-components=1 -C build/samtools

RUN cd build/samtools && \
    ./configure && \
    make install

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /usr/local/bin/samtools /usr/local/bin
