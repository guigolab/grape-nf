# Dockerfile for sambamba, compiled from source
#
FROM ubuntu:16.04 AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential \
                                               ca-certificates \
                                               curl \
                                               git \
                                               llvm \
                                               python3 \
                                               zlib1g-dev

ARG LDC_VERSION=1.11.0

RUN curl -fsSL "https://github.com/ldc-developers/ldc/releases/download/v${LDC_VERSION}/ldc2-${LDC_VERSION}-linux-x86_64.tar.xz" \
    | tar xJ && \
    mv "/ldc2-${LDC_VERSION}-linux-x86_64" "/ldc"

ENV \
    PATH="/ldc/bin:${PATH}" \
    LD_LIBRARY_PATH="/ldc/lib:/usr/lib:/lib:${LD_LIBRARY_PATH}" \
    LIBRARY_PATH="/ldc/lib:/usr/lib:/lib:${LD_LIBRARY_PATH}"

ARG REVISION=v0.7.0

RUN git clone --recursive https://github.com/biod/sambamba.git && \
    cd sambamba && git checkout $REVISION && \
    make

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /sambamba/bin/sambamba /usr/local/bin
