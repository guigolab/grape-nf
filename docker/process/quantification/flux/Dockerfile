# Dockerfile for the grape-nf FLUX-CAPACITOR quantification image
#
ARG SAMTOOLS_VER=1.3.1
ARG FLUX_VER=1.6.1

FROM grapenf/samtools:${SAMTOOLS_VER} AS samtools
FROM grapenf/flux-capacitor:${FLUX_VER}

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

ARG SAMTOOLS_VER=1.3.1
ENV _SAMTOOLS_VER ${SAMTOOLS_VER}

COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/
