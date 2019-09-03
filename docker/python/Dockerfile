# Dockerfile for the grape-nf python image
#
FROM       grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache py2-pip python2

RUN apk add --no-cache --virtual=.build-dependencies build-base \
                                             python-dev \
                                             zlib-dev \
                                             bzip2-dev \
                                             xz-dev && \
    pip --no-cache-dir install cython numpy pysam && \
    apk del .build-dependencies
