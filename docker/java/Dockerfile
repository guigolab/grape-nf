# Dockerfile for the grape-nf base image
#
FROM       frolvlad/alpine-oraclejre8

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache \
    bash \
    bzip2 \
    ca-certificates \
    coreutils \
    curl \
    findutils \
    libgomp \
    libstdc++ \
    pigz \
    procps \
    tar \
    xz-libs

COPY --from=ubuntu:16.04 /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/libz.so.1.2.8 /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21 /usr/glibc-compat/lib/

RUN \
# Fix ldd
    sed -i -e 's#/bin/bash#/bin/sh#' /usr/glibc-compat/bin/ldd && \
    sed -i -e 's#RTLDLIST=.*#RTLDLIST="/usr/glibc-compat/lib/ld-linux-x86-64.so.2"#' /usr/glibc-compat/bin/ldd && \
    \
# Isolate the glibc lib's from musl-libc lib's
    rm \
        /lib/ld-linux-x86-64.so.2 \
        /etc/ld.so.cache && \
    echo "/usr/glibc-compat/lib" > /usr/glibc-compat/etc/ld.so.conf && \
    /usr/glibc-compat/sbin/ldconfig -i && \
    rm -r \
        /var/cache/ldconfig
