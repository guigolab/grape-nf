# Dockerfile for sambamba
#
FROM grapenf/base

ARG SAMBAMBA_VER=0.7.0

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ENV _SAMBAMBA_VERSION $SAMBAMBA_VER

RUN curl -fsSL https://github.com/biod/sambamba/releases/download/v$_SAMBAMBA_VERSION/sambamba-${_SAMBAMBA_VERSION}-linux-static.gz \
    | pigz -dc > /usr/local/bin/sambamba && \
    chmod +x /usr/local/bin/sambamba