#!/usr/bin/awk -f
#
#  Copyright (c) 2013, Centre for Genomic Regulation (CRG)
#  Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
#
#  This file is part of the Blueprint RNAseq pipeline.
#
#  Ensure gff format correctness
#

$1!~/#/{
    for (i=1;i<=7;i++)
    {
	printf $i"\t";
    }
    printf $8;
    if(NF>8)
    {
	printf "\t"$9;
	for (i=10;i<=NF;i++)
	{
	    printf " "$i;
	}
    }
    print "";
}

$1~/#/{print}
