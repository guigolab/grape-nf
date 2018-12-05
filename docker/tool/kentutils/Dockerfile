# Dockerfile for Kent utils from UCSC
#
FROM ubuntu:16.04 AS builder

ENV _KENTUTILS_VERSION 308

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    unzip \
    libssl-dev \
    libpng-dev \
    libmysqlclient-dev

RUN mkdir build && \
    curl http://hgdownload.soe.ucsc.edu/admin/exe/userApps.archive/userApps.v$_KENTUTILS_VERSION.src.tgz \
    | tar xz -C build

RUN cd build/userApps && make libs && \
    cd kent/src && make destBin && \
    cd utils && make bedGraphToBigWig.userApps && \
    cd ../hg && make genePredToBed.userApp && \
    cd utils && make gtfToGenePred.userApp

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /usr/lib/x86_64-linux-gnu/libmysqlclient.so.20 /usr/glibc-compat/lib/
COPY --from=builder /lib/x86_64-linux-gnu/libssl.so.1.0.0 /usr/glibc-compat/lib/
COPY --from=builder /lib/x86_64-linux-gnu/libcrypto.so.1.0.0 /usr/glibc-compat/lib/
COPY --from=builder /build/userApps/bin/* /usr/local/bin/