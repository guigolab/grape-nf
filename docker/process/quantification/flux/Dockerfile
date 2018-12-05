# Dockerfile for the grape-nf FLUX-CAPACITOR quantification image
#
FROM grapenf/flux-capacitor:1.6.1

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=grapenf/samtools:1.3.1 /usr/local/bin/samtools /usr/local/bin/
COPY --from=grapenf/samtools:1.3.1 /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
