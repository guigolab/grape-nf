# Dockerfile for the grape-nf inferexp image
#
FROM grapenf/rseqc:2.6.4

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=grapenf/kentutils:308 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/kentutils:308 /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
