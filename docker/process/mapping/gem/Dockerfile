# Dockerfile for the grape-nf GEM mapping image
#
ARG GEMTOOLS_VER=1.7.1
ARG SAMTOOLS_VER=1.3.1

FROM grapenf/samtools:${SAMTOOLS_VER} as samtools

FROM grapenf/gemtools:${GEMTOOLS_VER}

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/
