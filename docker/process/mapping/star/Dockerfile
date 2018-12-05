# Dockerfile for the grape-nf STAR mapping step
#
FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=grapenf/star:2.4.0j /usr/local/bin/STAR /usr/local/bin/
COPY --from=grapenf/star:2.4.0j /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=grapenf/sambamba:0.6.8 /usr/local/bin/sambamba /usr/local/bin/
