# Dockerfile for the grape-nf contig/bigwig step
#
FROM grapenf/python

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=grapenf/bedtools:2.19.1 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/sambamba:0.6.8 /usr/local/bin/sambamba /usr/local/bin/
