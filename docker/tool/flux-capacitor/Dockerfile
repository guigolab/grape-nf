# Dockerfile for the Flux Capacitor
#
FROM grapenf/java

ARG FLUX_VER=1.6.1

ENV _FLUX_VERSION ${FLUX_VER}

ENV PATH $PATH:/flux-capacitor-1.6.1/bin

RUN curl http://artifactory.sammeth.net/artifactory/barna/barna/barna.capacitor/$_FLUX_VERSION/flux-capacitor-$_FLUX_VERSION.tgz | tar xzf -
