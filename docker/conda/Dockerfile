# Dockerfile for the grape-nf base image using conda
#
FROM frolvlad/alpine-miniconda3

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache bash \
                       coreutils \
                       gawk \
                       pigz

ADD grape-env.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
ENV PATH /opt/conda/envs/grape/bin:$PATH

RUN conda clean -ay