# Dockerfile for bamstats
#
FROM       grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ENV _BAMSTATS_VERSION 0.3.2

RUN curl -fsSL https://github.com/guigolab/bamstats/releases/download/v$_BAMSTATS_VERSION/bamstats-v$_BAMSTATS_VERSION-linux-amd64.tar.bz2 \
    | tar xj --strip-components 3 -C /usr/local/bin