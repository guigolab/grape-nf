# Dockerfile for the grape-nf contig/bigwig step
#
ARG BEDTOOLS_VER=2.19.1
ARG SAMBAMBA_VER=0.7.1

FROM grapenf/bedtools:${BEDTOOLS_VER} as bedtools
FROM grapenf/sambamba:${SAMBAMBA_VER} as sambamba

FROM grapenf/python

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=bedtools /usr/local/bin/* /usr/local/bin/
COPY --from=sambamba /usr/local/bin/sambamba /usr/local/bin/
