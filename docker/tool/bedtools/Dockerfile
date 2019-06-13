# Dockerfile for BEDtools
#
FROM ubuntu:16.04 AS builder

ARG BEDTOOLS_VER=2.19.1

ENV _BEDTOOLS_VERSION ${BEDTOOLS_VER}

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    curl \
    python \
    ca-certificates

RUN mkdir -p build/bedtools && \
    curl -L https://api.github.com/repos/arq5x/bedtools2/tarball/v$_BEDTOOLS_VERSION \
    | tar xzf - --strip-components=1 -C build/bedtools

RUN cd build/bedtools && make

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /build/bedtools/bin/* /usr/local/bin/