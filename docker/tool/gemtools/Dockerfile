# Dockerfile for GEMtools
#
FROM grapenf/python

ENV _GEMTOOLS_VERSION 1.7.1

ENV PATH /usr/lib/python2.7/site-packages/gem/gembinaries/:${PATH}

RUN apk add --no-cache --virtual=.build-dependencies \
        build-base \
        python-dev \
        zlib-dev \
        bzip2-dev \
        xz-dev && \
    pip --no-cache-dir install argparse cython==0.18

RUN mkdir -p build/gemtools && \
    curl -L https://api.github.com/repos/gemtools/gemtools/tarball/v${_GEMTOOLS_VERSION} \
    | tar xz --strip-components=1 -C build/gemtools && \
    sed -i '32 s/GENERAL_FLAGS=/GENERAL_FLAGS=-std=gnu89 /g' build/gemtools/GEMTools/Makefile.mk && \
    sed -i 's/http/https/g' build/gemtools/distribute_setup.py && \
    sed -i "s/\['z', 'bz2', 'gemtools'\]/\['gemtools', 'z', 'bz2'\]/" build/gemtools/setup.py && \
    sed -i 's/en_US.UTF-8/C/g' build/gemtools/python/gem/reports.py && \
    cd build/gemtools && python setup.py install

RUN apk del .build-dependencies && \
    # pip uninstall -y numpy && \
    rm -rf /build