# Dockerfile for the grape-nf RSEM quantification image
#
ARG SAMBAMBA_VER=0.7.0
ARG RSEM_VER=1.2.21

FROM grapenf/sambamba:${SAMBAMBA_VER} AS sambamba
FROM grapenf/rsem:${RSEM_VER} AS rsem

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache perl R

ARG SAMBAMBA_VER=0.7.0
ENV _SAMBAMBA_VER $SAMBAMBA_VER

COPY --from=rsem /usr/local/bin/* /usr/local/bin/
COPY --from=sambamba /usr/local/bin/sambamba /usr/local/bin/
