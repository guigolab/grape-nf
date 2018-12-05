# Dockerfile for samtools
#
FROM ubuntu:16.04 AS builder

ENV _SAMTOOLS_VERSION 1.3.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libncurses-dev \
    curl \
    ca-certificates

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
COPY --from=builder /lib/x86_64-linux-gnu/libncurses.so.5 /lib/x86_64-linux-gnu/libtinfo.so.5 /usr/glibc-compat/lib/