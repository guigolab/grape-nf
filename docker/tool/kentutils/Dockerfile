# Dockerfile for Kent utils from UCSC
#
ARG ALPINE_VER=3.11
FROM alpine:${ALPINE_VER} AS builder

ARG KENTUTILS_VER=308

ENV _KENTUTILS_VERSION ${KENTUTILS_VER}

RUN apk update && apk --no-cache add \
    build-base \
    curl \
    zlib-dev \
    openssl-dev \
    libpng-dev \
    mariadb-dev

RUN mkdir build && \
    curl http://hgdownload.soe.ucsc.edu/admin/exe/userApps.archive/userApps.v$_KENTUTILS_VERSION.src.tgz \
    | tar xz -C build

RUN cd build/userApps && \
    sed -i '1i#include <stdint.h>' samtabix/knetfile.h && \
    sed -i '11istruct _IO_FILE { char _x; };' kent/src/lib/fof.c && make libs && \
    cd kent/src && make destBin && \
    cd utils && make bedGraphToBigWig.userApps && \
    cd ../hg && make genePredToBed.userApp && \
    cd utils && make gtfToGenePred.userApp

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk update && \
    apk --no-cache add mariadb-connector-c

COPY --from=builder /build/userApps/bin/* /usr/local/bin/
