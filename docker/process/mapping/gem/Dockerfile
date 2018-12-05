# Dockerfile for the grape-nf GEM mapping image
#
FROM grapenf/gemtools:1.7.1

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=grapenf/samtools:1.3.1 /usr/local/bin/samtools /usr/local/bin/
COPY --from=grapenf/samtools:1.3.1 /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
