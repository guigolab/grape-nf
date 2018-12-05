# Dockerfile for GEMtools
#
FROM ubuntu:16.04 AS builder

ENV _GEMTOOLS_VERSION 1.7.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        curl \
        python \
        python-pip \
        python-dev \
        ca-certificates && \
    pip --no-cache-dir install cython==0.18 

RUN mkdir -p build/gemtools && \
    curl -L https://api.github.com/repos/gemtools/gemtools/tarball/v${_GEMTOOLS_VERSION} \
    | tar xz --strip-components=1 -C build/gemtools && \
    sed -i '32 s/GENERAL_FLAGS=/GENERAL_FLAGS=-std=gnu89 /g' build/gemtools/GEMTools/Makefile.mk && \
    sed -i 's/http/https/g' build/gemtools/distribute_setup.py && \
    sed -i "s/\['z', 'bz2', 'gemtools'\]/\['gemtools', 'z', 'bz2'\]/" build/gemtools/setup.py && \
    sed -i 's/en_US.UTF-8/C/g' build/gemtools/python/gem/reports.py && \
    cd build/gemtools &&\
    python distribute_setup.py && \
    make dist

FROM grapenf/python

ENV PATH /opt/gemtools/bin:/opt/gentools/lib/python2.7/site-packages/gem/gembinaries/:${PATH}

COPY --from=builder /build/gemtools/dist/gemtools*i3 /opt/gemtools
