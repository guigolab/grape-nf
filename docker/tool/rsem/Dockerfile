# Dockerfile for RSEM
#
FROM ubuntu:16.04 AS builder

ENV _RSEM_VERSION 1.2.21

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libncurses-dev \
    curl \
    ca-certificates

RUN mkdir -p build/rsem && curl -L http://deweylab.biostat.wisc.edu/rsem/src/rsem-$_RSEM_VERSION.tar.gz | tar xzf - --strip-components=1 -C build/rsem

RUN cd build/rsem && make

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache perl R

COPY --from=builder /build/rsem/rsem-extract-reference-transcripts /usr/local/bin
COPY --from=builder /build/rsem/rsem-synthesis-reference-transcripts /usr/local/bin
COPY --from=builder /build/rsem/rsem-preref /usr/local/bin
COPY --from=builder /build/rsem/rsem-parse-alignments /usr/local/bin
COPY --from=builder /build/rsem/rsem-build-read-index /usr/local/bin
COPY --from=builder /build/rsem/rsem-run-em /usr/local/bin
COPY --from=builder /build/rsem/rsem-tbam2gbam /usr/local/bin
COPY --from=builder /build/rsem/rsem-run-gibbs /usr/local/bin
COPY --from=builder /build/rsem/rsem-calculate-credibility-intervals /usr/local/bin
COPY --from=builder /build/rsem/rsem-simulate-reads /usr/local/bin
COPY --from=builder /build/rsem/rsem-bam2wig /usr/local/bin
COPY --from=builder /build/rsem/rsem-get-unique /usr/local/bin
COPY --from=builder /build/rsem/rsem-bam2readdepth /usr/local/bin
COPY --from=builder /build/rsem/rsem-sam-validator /usr/local/bin
COPY --from=builder /build/rsem/rsem-scan-for-paired-end-reads /usr/local/bin
COPY --from=builder /build/rsem/rsem_perl_utils.pm /usr/local/bin
COPY --from=builder /build/rsem/convert-sam-for-rsem /usr/local/bin
COPY --from=builder /build/rsem/extract-transcript-to-gene-map-from-trinity /usr/local/bin
COPY --from=builder /build/rsem/rsem-calculate-expression /usr/local/bin
COPY --from=builder /build/rsem/rsem-control-fdr /usr/local/bin
COPY --from=builder /build/rsem/rsem-gen-transcript-plots /usr/local/bin
COPY --from=builder /build/rsem/rsem-generate-data-matrix /usr/local/bin
COPY --from=builder /build/rsem/rsem-generate-ngvector /usr/local/bin
COPY --from=builder /build/rsem/rsem-plot-model /usr/local/bin
COPY --from=builder /build/rsem/rsem-plot-transcript-wiggles /usr/local/bin
COPY --from=builder /build/rsem/rsem-prepare-reference /usr/local/bin
COPY --from=builder /build/rsem/rsem-run-ebseq /usr/local/bin
