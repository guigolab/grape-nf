#!/bin/bash 
set -e
set -o pipefail

RED="\033[1;31m"
GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
NORMAL="\033[0m"

[ ! -z "$1" ] && config="-c $1"
[ ! -z "$2" ] && chunks="--chunk-size $2"

echo -e "==$YELLOW Running pipeline$NORMAL"
nextflow run grape-pipeline.nf --index test-index.txt --genome data/genome.fa --annotation data/annotation.gtf -process.executor=local ${config} ${chunks} | tee test.log
echo -e "==$YELLOW Compare results$NORMAL"
dbFile=$(grep "db ->" test.log | cut -d " " -f5) || (echo -e "==$RED [ERROR] pipeline db file path not found!$NORMAL" >>/dev/stderr; exit 1)
cat <(cat $dbFile | cut -f1 | xargs md5sum) data/refout.md5 | sort | awk '{sums[$1]++; n=split($2,a,"/");files[$1]=a[n]}END{exitStatus=0;for(sum in sums){printf("%-32s\t%2s\n", files[sum], sums[sum]!=2 ? "\033[1;31mKO\033[0m" : "\033[1;32mOK\033[0m"); if(sums[sum]!=2){exitStatus=1}}; exit exitStatus}'
echo -e "==$YELLOW DONE$NORMAL"
