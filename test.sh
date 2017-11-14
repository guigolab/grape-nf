#!/bin/bash 
set -e
set -o pipefail

usage() {
	echo "Usage:"
	echo "  $0 <profile>"
	exit 1
}

mkcd() {
	local DIR=$1
	mkdir -p $DIR && cd $DIR
}

RED="\033[1;31m"
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
NORMAL="\033[0m"

[[ $1 == "help" ]] && usage

PROFILE=${PROFILE-"testGemFlux"} 
CHECKDIR=${CHECKDIR-"checksum"}

RUN_OPTS=${RUN_OPTS-"-process.errorStrategy=terminate"}

echo -e "==$YELLOW Running pipeline with profile -> $BLUE${PROFILE}$NORMAL"
[ ! -x nextflow ] && (curl -fsSL get.nextflow.io | bash && chmod +x nextflow) || true
./nextflow -c test-profiles.config run . -profile $PROFILE ${RUN_OPTS}
echo -e "==$YELLOW Compare results$NORMAL"
mkcd $CHECKDIR
cut -f 3 ../pipeline.db | xargs -I{} ln -fs {}
md5sum -c ../data/$PROFILE.md5
echo -e "==$YELLOW DONE$NORMAL"
