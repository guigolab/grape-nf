# Dockerfile for the grape-nf contig/bigwig step
#
ARG KENTUTILS_VER=308
ARG SAMBAMBA_VER=0.7.0
ARG BEDTOOLS_VER=2.19.1

FROM grapenf/bedtools:${BEDTOOLS_VER} as bedtools
FROM grapenf/kentutils:${KENTUTILS_VER} as kentutils
FROM grapenf/sambamba:${SAMBAMBA_VER} as sambamba

FROM grapenf/python

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=bedtools /usr/local/bin/* /usr/local/bin/
COPY --from=kentutils /usr/local/bin/* /usr/local/bin/
COPY --from=kentutils /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=sambamba /usr/local/bin/sambamba /usr/local/bin/