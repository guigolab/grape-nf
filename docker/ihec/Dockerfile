# Dockerfile for the grape-nf IHEC image
#
ARG STAR_VER=2.4.0j
ARG KENTUTILS_VER=308
ARG SAMBAMBA_VER=0.7.1
ARG BAMSTATS_VER=0.3.4
ARG RSEM_VER=1.2.21
ARG RSEQC_VER=2.6.4
ARG BEDTOOLS_VER=2.19.1
ARG SAMTOOLS_VER=1.3.1

FROM grapenf/star:${STAR_VER} as star
FROM grapenf/bedtools:${BEDTOOLS_VER} as bedtools
FROM grapenf/kentutils:${KENTUTILS_VER} as kentutils
FROM grapenf/sambamba:${SAMBAMBA_VER} as sambamba
FROM grapenf/samtools:${SAMTOOLS_VER} as samtools
FROM grapenf/bamstats:${BAMSTATS_VER} as bamstats
FROM grapenf/rsem:${RSEM_VER} as rsem

FROM grapenf/rseqc:${RSEQC_VER}

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache perl R mariadb-connector-c

COPY --from=star /usr/local/bin/STAR /usr/local/bin/
COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/
COPY --from=bedtools /usr/local/bin/* /usr/local/bin/
COPY --from=kentutils /usr/local/bin/* /usr/local/bin/
COPY --from=rsem /usr/local/bin/* /usr/local/bin/
COPY --from=sambamba /usr/local/bin/sambamba /usr/local/bin/
COPY --from=bamstats /usr/local/bin/bamstats /usr/local/bin/
