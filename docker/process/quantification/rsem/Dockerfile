# Dockerfile for the grape-nf RSEM quantification image
#
FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

RUN apk add --no-cache perl R

COPY --from=grapenf/rsem:1.2.21 /usr/local/bin/* /usr/local/bin/
COPY --from=grapenf/sambamba:0.6.8 /usr/local/bin/sambamba /usr/local/bin/
