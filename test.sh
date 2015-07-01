#!/bin/bash 
set -e
set -o pipefail

usage() {
	echo "Usage:"
	echo "  $0 <profile> <nextflow_args>"
	exit 1
}

RED="\033[1;31m"
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
NORMAL="\033[0m"

[ $# -eq 0 ] && usage

profile=${1} && shift
[ -z $profile ] && profile="gemflux"

log_file="test-$profile.log"

nxf_opts="-process.executor=local -process.errorStrategy=terminate ${@}"

echo -e "==$YELLOW Running pipeline with $BLUE<${profile}>$YELLOW profile$NORMAL"
nextflow run grape-pipeline.nf -profile $profile --index test-index.txt --genome data/genome.fa --annotation data/annotation.gtf ${nxf_opts} > $log_file
echo -e "==$YELLOW Compare results$NORMAL"
dbFile=$(grep "db ->" $log_file | cut -d " " -f5) || (echo -e "==$RED [ERROR] pipeline db file path not found!$NORMAL" >>/dev/stderr; exit 1)
cat <(cat $dbFile | cut -f3 | xargs md5sum) data/${profile}.md5 | sort | awk '{sums[$1]++; n=split($2,a,"/");files[$1]=a[n]}END{exitStatus=0;for(sum in sums){printf("%-32s\t%2s\n", files[sum], sums[sum]!=2 ? "\033[1;31mKO\033[0m" : "\033[1;32mOK\033[0m"); if(sums[sum]!=2){exitStatus=1}}; exit exitStatus}'
echo -e "==$YELLOW DONE$NORMAL"
