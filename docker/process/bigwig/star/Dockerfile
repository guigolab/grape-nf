# Dockerfile for the grape-nf STAR mapping step
#
ARG STAR_VER=2.4.0j
ARG KENTUTILS_VER=308
ARG SAMBAMBA_VER=0.7.0

FROM grapenf/star:${STAR_VER} as star
FROM grapenf/kentutils:${KENTUTILS_VER} as kentutils
FROM grapenf/sambamba:${SAMBAMBA_VER} as sambamba

FROM grapenf/base

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=star /usr/local/bin/STAR /usr/local/bin/
COPY --from=star /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=kentutils /usr/local/bin/* /usr/local/bin/
COPY --from=kentutils /usr/glibc-compat/lib/* /usr/glibc-compat/lib/
COPY --from=sambamba /usr/local/bin/sambamba /usr/local/bin/