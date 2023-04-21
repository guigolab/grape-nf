# Dockerfile for samtools
#
FROM frolvlad/alpine-miniconda3

LABEL author.name="Paolo Marangio"
LABEL author.email="marangiopaolo@gmail.com"

RUN apk update && apk add --no-cache bash \
                                    coreutils \
                                    gawk \
                                    pigz



ADD kallisto-env.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
ENV PATH /opt/conda/envs/kallisto/bin:$PATH

RUN conda clean -ay
