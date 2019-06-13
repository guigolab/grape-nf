# Dockerfile for bamstats
#
FROM grapenf/base

ARG BAMSTATS_VER=0.3.2

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ENV _BAMSTATS_VERSION ${BAMSTATS_VER}

RUN curl -fsSL https://github.com/guigolab/bamstats/releases/download/v$_BAMSTATS_VERSION/bamstats-v$_BAMSTATS_VERSION-linux-amd64.tar.bz2 \
    | tar xj --strip-components 3 -C /usr/local/bin