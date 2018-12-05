# Dockerfile for STAR
#
FROM ubuntu:16.04 AS builder

ENV _STAR_VERSION 2.4.0j

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    curl \
    ca-certificates \
    vim-common  

RUN mkdir -p build/STAR && curl -L https://api.github.com/repos/alexdobin/STAR/tarball/STAR_$_STAR_VERSION | tar xz --strip-components=1 -C build/STAR 

RUN cd build/STAR/source && make STAR

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /usr/lib/x86_64-linux-gnu/libgomp.so.1 /usr/glibc-compat/lib/
COPY --from=builder /build/STAR/source/STAR /usr/local/bin/

