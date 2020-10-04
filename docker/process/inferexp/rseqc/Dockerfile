# Dockerfile for the grape-nf inferexp image
#
ARG RSEQC_VER=2.6.4
ARG KENTUTILS_VER=308

FROM grapenf/kentutils:${KENTUTILS_VER} as kentutils

FROM grapenf/rseqc:${RSEQC_VER}

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache mariadb-connector-c

COPY --from=kentutils /usr/local/bin/* /usr/local/bin/
