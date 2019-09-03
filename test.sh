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

canonicalize_path() {
    if [ -d "$1" ]; then
        _canonicalize_dir_path "$1"
    else
        _canonicalize_file_path "$1"
    fi
}   

_canonicalize_dir_path() {
    (cd "$1" 2>/dev/null && pwd -P) 
}           

_canonicalize_file_path() {
    local dir file
    dir=$(dirname -- "$1")
    file=$(basename -- "$1")
    (cd "$dir" 2>/dev/null && printf '%s/%s\n' "$(pwd -P)" "$file")
}

nxf_setup() {
	if ! nextflow &>/dev/null; then
		export PATH=$PWD:$PATH 
		if [ ! -x nextflow ]; then
			curl -fsSL get.nextflow.io | bash && chmod +x nextflow
		fi
	fi
}

normalize_bam() {
	docker run --rm -v $1:$1 -w $PWD grapenf/samtools:1.3.1 bash -c "samtools view -h $2 | grep -v '@CO\|@PG' | samtools view -Sb - > $(basename $2) 2>/dev/null"
}

get_file() {
	baseDir=$1
	filename=$2
	case "$filename" in
    	*.bam) normalize_bam "$baseDir" "$filename" ;;
    	*)         ln -fs "$filename"
	esac
}

export -f get_file normalize_bam

RED="\033[1;31m"
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
NORMAL="\033[0m"

[[ $1 == "help" ]] && usage

BASE_DIR=$(canonicalize_path $(dirname $0))

PROFILE=${PROFILE-"starrsem"} 
CHECKDIR=${CHECKDIR-"checksum"}

RUN_OPTS=${RUN_OPTS-"-process.errorStrategy=terminate"}

WITH_DOCKER=${WITH_DOCKER-"1"}
[[ $WITH_DOCKER ]] && RUN_OPTS=" -with-docker"

echo -e "==$YELLOW Running pipeline with profile -> $BLUE${PROFILE}$NORMAL"
nxf_setup && nextflow run ${BASE_DIR} -profile $PROFILE ${RUN_OPTS} "$@"
echo -e "==$YELLOW Compare results$NORMAL"
mkcd $CHECKDIR
cut -f 3 ../pipeline.db | xargs -I{} bash -c "get_file $BASE_DIR {}"
md5sum -c ${BASE_DIR}/data/$PROFILE.md5
echo -e "==$YELLOW Clean up$NORMAL"
cd .. && rm -rf $CHECKDIR
echo -e "==$YELLOW DONE$NORMAL"
