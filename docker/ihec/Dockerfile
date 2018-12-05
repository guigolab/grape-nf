# Dockerfile for the grape-nf IHEC image
#
FROM grapenf/rseqc:2.6.4

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache perl R

COPY --from=grapenf/star:2.4.0j /usr/local/bin/STAR /usr/local/bin/
COPY --from=grapenf/star:2.4.0j /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=grapenf/samtools:1.3.1 /usr/local/bin/samtools /usr/local/bin/
COPY --from=grapenf/samtools:1.3.1 /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=grapenf/bedtools:2.19.1 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/kentutils:308 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/kentutils:308 /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=grapenf/rsem:1.2.21 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/sambamba:0.6.8 /usr/local/bin/sambamba /usr/local/bin/
COPY --from=grapenf/bamstats:0.3.2 /usr/local/bin/bamstats /usr/local/bin/
