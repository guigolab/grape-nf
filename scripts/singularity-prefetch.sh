#!/bin/bash
set -e
set -u

usage() {
    cat << USAGE
Usage: $0 [OPTION]...
Pull the IHEC singularity image into 'NXF_SINGULARITY_CACHEDIR' ('work/singularity' if the environment variable is not defined).
By default, it does not overwrite the destination image files if it already exists.

Optional arguments:
-f, --force    Force overwriting destination image if it exists.
-h, --help     Print usage information and exit.
USAGE
    exit 1
}

# Docker imaged to be pulled
imgUrls=(
    grapenf/ihec:latest
)

# use NXF_SINGULARITY_CACHEDIR if defined or 'work/singularity'
outdir=${NXF_SINGULARITY_CACHEDIR-work/singularity}
mkdir -p $outdir

# check arguments
FORCE=
case ${1-} in 
    -f | --force)
        FORCE=-F
        shift;;
    -h | --help)
        usage
        shift;;
esac

for img in ${imgUrls}; do
    name=$outdir/${img//[:\/]/-}.img
    singularity pull $FORCE --name ${name} docker://$img 2>/dev/null \
        || echo "Skipping $name. Image already exists."
done