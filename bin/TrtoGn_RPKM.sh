#!/bin/bash
#
#  Copyright (c) 2013, Centre for Genomic Regulation (CRG)
#  Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
#
#  This file is part of the Blueprint RNAseq pipeline.
#

# takes as input a gtf file of transcripts produced by the flux capacitor
# and makes as output a gene gff file with an rpkm value associated to each gene
# and computed as the sum of the rpkm of the transcripts belonging to this gene
# as well as the number of reads supporting the gene computed as the sum of the
# reads supporting each transcript of this gene.
# Be careful: makes the following assumptions:
##############################################
# 1. the annotation file has gene_id gene as first key/value pair, and transcript_id transcript as second key/value pair
# 2. the flux output file has transcript_id transcript as first key/value pair

# Programs
##########
prog=`basename $0`
gff2gff="gff2gff.awk"

# Usage
#######
function usage {    
    echo "" >&2
    echo "Produce a gene gff file with an rpkm value and a number of read value for each gene, corresponding" >&2
    echo "to the sum of the rpkm and of the number of reads of the transcripts belonging to this gene." >&2
    echo "" >&2
    echo "Usage: $prog -a ANNOTATION -i TRANSCRIPTS_GTF -o OUTPUT_FOLDER" >&2
    echo "" >&2
    printf "\t-a|--annnotation\t an annotation gtf file\n" >&2
    printf "\t-i|--input\t\t a flux output file with rpkm value\n" >&2
    echo "" >&2
    echo "Options:"
    echo "" >&2
    printf "\t-o|--output\t\t the output folder and filename. Default: ./<input_basename>_gene_with_rpkm.gff\n"
    echo "" >&2
    exit 1
}

# Execute getopt
ARGS=`getopt -o "i:a:o:h" -l "input:,annotation:,output:,help" \
      -n "$prog" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

# A little magic
eval set -- "$ARGS"

output="."

# Now go through all the options
while true;
do
  case "$1" in
    -i|--input)
      if [ -n "$2" ];
      then
        tr=$2
      fi
      shift 2;;

    -a|--annotation)
      if [ -n "$2" ];
      then
        annot=$2
      fi
      shift 2;;

    -o|--output)
      if [ -n "$2" ];
      then
        output=$2
      fi
      shift 2;;
    
    -h|--help)
      usage
      shift;;

    --)
      shift
      break;;
  esac
done

if [ ! $tr ] || [ ! $annot ];then
    usage
    exit 1
fi

annotbase=`basename $annot .gz`
trbase=`basename $tr`
[ output == '.' ] && output=${trbase%.gtf}\_gene_with_rpkm.gff
withtrlist=$(dirname $output)/${trbase}.${annotbase%.gtf}.withtrlist.gff

cat=cat
if [[ $annot =~ \.gz$ ]]; then
  cat=zcat
fi

echo "I am making the file of genes with associated transcripts from the annotation" >&2
$cat $annot | awk 'BEGIN{OFS=FS="\t"}$3=="transcript"{
	split($9, a, "; "); for(i=1;i<=length(a);i++){split(a[i], b, " "); gsub(/"/, "", b[2]); dict[b[1]]=b[2]}
	gene_id=dict["gene_id"]; transcript_id=dict["transcript_id"]; 
	trlist[gene_id]=(trlist[gene_id])(transcript_id)(",");}
$3=="gene"{
	split($9, a, "; "); for(i=1;i<=length(a);i++){split(a[i], b, " "); gsub(/"/, "", b[2]); dict[b[1]]=b[2]}; 
	gene_id=dict["gene_id"]; transcript_id=dict["transcript_id"]; 
	line[gene_id]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \""gene_id"\""}
END{for(g in line){print line[g]"; transcript_ids \""trlist[g]"\";"}}' | $gff2gff > $withtrlist

echo "I am making the gene file with rpkm and number of reads" >&2
awk -v fileRef=$tr 'BEGIN{while (getline < fileRef >0){
  k=9; while(k<=(NF-1)){
      split($10,a,"\""); 
      if($k=="RPKM"){split($(k+1),b,";"); rpkm[a[2]]=b[1];} 
      if($k=="reads"){split($(k+1),b,";"); reads[a[2]]=b[1];} k+=2}
    }
  } {
    split($12,a,"\""); split(a[2],b,","); s1=0; k=1; while(b[k]!=""){
      s1+=rpkm[b[k]]; k++
    } s2=0; k=1; while(b[k]!=""){
      s2+=reads[b[k]]; k++
    } print $0, "RPKM", s1";", "reads", s2";"
  }' $withtrlist | $gff2gff | sort -k1,1 -k4,4n > $output

echo "I am removing unuseful files" >&2
rm $withtrlist

echo "I am done" >&2


