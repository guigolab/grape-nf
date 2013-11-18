#!/bin/bash
#
#Forcing bash shell
#
#$ -S /bin/bash
#
#$ -V
#
#$ -cwd
# -M emiliopalumbo@gmail.com
# -m b
#
#$ -pe smp 8
#$ -q rg-el6,long
#$ -l virtual_free=64G,h_rt=240:00:00
#
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err
#

#set -eo pipefail

function usage {
    echo ""
    echo "### Blueprint RNAseq pipeline ###"
    echo "Run the RNAseq pipeline on one sample."
    echo ""
    echo "Usage: $0 -i FASTQ_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]..."
    echo ""
    printf "  -i|--input\t\tinput file\n"
    printf "  -g|--genome\t\treference genome file\n"
    printf "  -a|--annotation\treference gene annotation file\n"
    echo ""
    echo "Options:"
    printf "  -m|--mismatches\tMax number of mismatches. Default \"4\"\n"
    printf "  -n|--hits\t\tMax number of hits. Default \"10\"\n"
    printf "  -q|--quality-offset\tThe quality offset of the fastq files. Default: \"33\".\n"    
    printf "  -r|--max-read-length\tThe maximum read length (used to compute the transcriptomes). Default: \"150\".\n"    
    printf "  -s|--read-strand\tdirectionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default \"NONE\".\n"
    printf "  -l|--loglevel\t\tLog level (error, warn, info, debug). Default \"info\".\n"
    printf "  -t|--tmp-dir\t\tSpecify local temporary folder to copy files when running on distributed file systems. Default: no tmp folder.\n"
    printf "  -h|--help\t\tShow this message and exit.\n"
    printf "  --dry\t\t\tTest the pipeline. Writes the command to the standard output.\n"
    exit 0
}

function printHeader {
    local string=$1
    echo "`date` *** $string ***"
}

function log {
    local string=$1
    local label=$2
    if [[ ! $ECHO ]];then
        if [[ $label != "" ]];then
            printf "[$label] $string"
        else
            printf "$string"
        fi
    fi
}

function getAbsPath {
    local path=$1
    echo "`readlink -en $path`"
}

function run {
    local command=($1)
    if [[ $2 ]];then
        ${2}${command[@]}
    else
        eval ${command[@]}
    fi
}

function copyToTmp {    
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        local name=`basename $i`
        if [ ! -e $tmpdir/$name ];then
            log "Copying $i to $TMPDIR..." $step
            run "cp $i $tmpdir" "$ECHO"
            log "done\n"
        fi
    done
}

## Parsing arguments
#

# Execute getopt
ARGS=`getopt -o "i:g:a:m:n:s:t:l:q:r:h" -l "input:,genome:,annotation:,mismatches:,hits:,read-strand:,threads:,loglevel:,quality:,max-read-length:,tmp-dir:,dry-run,help" \
      -n "run.pipeline.sh" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

# A little magic
eval set -- "$ARGS"

# Setting defaults

# pipeline
mism=4
hits=10
readStrand="NONE"
qualityOffset="33"
maxReadLength="150"

# general
loglevel="info"
threads="1"
tmpdir=${TMPDIR-"-"}

while true;
do
  case "$1" in
    -i|--input)
      if [ -n "$2" ];
      then
        input=$(getAbsPath $2)
      fi
      shift 2;;

    -g|--genome)
      if [ -n "$2" ];
      then
        genome=$(getAbsPath $2)
      fi
      shift 2;;

    -a|--annotation)
      if [ -n "$2" ];
      then
        annotation=$(getAbsPath $2)
      fi
      shift 2;;
 
    -m|--mismatches)
       if [ -n "$2" ];
       then
         mism=$2
       fi
       shift 2;;
    
    -n|--hits)
       if [ -n "$2" ];
       then
         hits=$2
       fi
       shift 2;;
    
    -s|--read-strand)
      if [ -n $2 ];
      then
        readStrand=$2
      fi
      shift 2;;
     
    -q|--quality)
      if [ -n $2 ];
      then
        qualityOffset=$2
      fi
      shift 2;;
     
    -r|--max-read-length)
      if [ -n $2 ];
      then
        maxReadLength=$2
      fi
      shift 2;;

    -l|--loglevel)
      if [ -n $2 ];
      then
        loglevel=$2
      fi
      shift 2;;
 
    -t|--threads)
      if [ -n $2 ];
      then
        threads=$2
      fi
      shift 2;;
 
    --tmp-dir)
      if [ -n $2 ];
      then
        tmpdir=$(getAbsPath $2)
      fi
      shift 2;;
    
    --dry-run)
      ECHO="echo "
      shift;;
    
    -h|--help)
      usage
      shift;;
    
    --)
      shift
      break;;
  esac
done

# Setting up environment

BASEDIR=`dirname ${SGE_O_WORKDIR-$PWD}`
BINDIR="$BASEDIR/bin"
export PATH=$BASEDIR/gemtools-1.6.2-i3/bin:$BASEDIR/flux-capacitor-1.2.4/bin:$BINDIR:$PATH

## Setting variables and input files
##
if [[ $input == "" ]];then
    log "Please specify the input file\n" "ERROR" >&2
    exit -1
fi

if [[ $genome == "" ]];then
    log "Please specify the genome file\n" "ERROR" >&2
    exit -1
fi

if [[ $annotation == "" ]];then
    log "Please specify the annotation file\n" "ERROR" >&2
    exit -1
fi

basename=$(basename $input)
sample=${basename%[_-\.]1*}

index="$genome.gem"

annName=`basename $annotation`

hthreads=$((threads/2))
if [[ $hthreads == 0 ]];then
    hthreads=1
fi

## Binaries
#
gem2sam="gem-2-sam"
samtools="samtools"
addXS="sam2cufflinks.sh"
trToGn="TrtoGn_RPKM.sh"
trToEx="TrtoEx_RPKM.sh"
bamToContigs="bamToContigs.sh"
gt_quality="gt.quality"
gt_filter="gt.filter"
gt_stats="gt.stats"
pigz="pigz"
BAMFLAG="bamflag"
makecontig="contigsNew.py"

## Print pipeline configuration

header="Pipeline configuration for $sample"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n" "Input file:" "$input"
printf "  %-34s %s\n" "Reference genome file:" "$genome"
printf "  %-34s %s\n" "Reference gene annotation file:" "$annotation"
printf "  %-34s %s\n" "Max number of allowed mismatches:" "$mism"
printf "  %-34s %s\n" "Max number of hits:" "$hits"
printf "  %-34s %s\n" "Quality offset:" "$qualityOffset"
printf "  %-34s %s\n" "Max read length:" "$maxReadLength"
printf "  %-34s %s\n" "Strandedness:" "$readStrand"
printf "  %-34s %s\n" "Number of threads:" "$threads"
printf "  %-34s %s\n" "Temporary folder:" "$tmpdir"
printf "  %-34s %s\n" "Loglevel:" "$loglevel"
echo ""

## START
#

printHeader "Starting Blueprint pipeline for $sample"
pipelineStart=$(date +%s)

## Mapping
#
if [[ `basename $input` =~ fastq ]];then
    if [ ! -e $sample.map.gz ];then
        step="MAP"
        startTime=$(date +%s)
        printHeader "Executing mapping step"
    
        ## Copy needed files to TMPDIR
        copyToTmp "index,annotation,t-index,keys"
    
        log "Running gemtools rna pipeline on ${sample}" $step
        run "gemtools --loglevel $loglevel rna-pipeline -f $input -q 33 -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -t $threads --no-bam" "$ECHO"
        #gemtools --loglevel $loglevel rna-pipeline -f $TMPDIR/$basename -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -t $threads -o $TMPDIR --no-bam
        #gemtools --loglevel $loglevel rna-pipeline -f $TMPDIR/$basename -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -m 150 -t $threads -o $TMPDIR --no-sam
    
        if [ -f $TMPDIR/${sample}.map.gz ]; then
            log "Computing md5sum for map file..." $step
            run "md5sum $TMPDIR/$sample.map.gz > $TMPDIR/$sample.map.gz.md5" "$ECHO"
            run "cp $TMPDIR/$sample.map.gz.md5 ." "$ECHO"
            log "done\n"
            log "Copying map file..." $step
            run "cp $TMPDIR/${sample}.map.gz ." "$ECHO"
            log "done\n"
        #else
        #    log "Error producing map file" "ERROR" >&2
        #    exit -1
        fi
        endTime=$(date +%s)
        printHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Map file already present...skipping mapping step"
    fi
    
  
    ## Filtering the map file
    ##
    filteredGem=${sample}_mism_${mism}_mmaps.map.gz
    
    if [ ! -e $filteredGem ];then
        step="FILTER"
        startTime=$(date +%s)
        printHeader "Executing filtering step"
    
        log "Filtering map file..." $step
        run "$gt_quality -i $sample.map.gz -t $threads | $gt_filter --max-levenshtein-error $mism -t $threads | $gt_filter --max-matches 10 -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
        log "done\n" $step
        if [ -f $filteredGem ]; then
            log "Computing md5sum for filtered file..." $step
            run "md5sum $filteredGem > $filteredGem.md5" "$ECHO"
            log "done\n"
        # else
        #    log "Error producing filtered map file" "ERROR" >&2
        #    exit -1
        fi
        endTime=$(date +%s)
        printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Filtered map file is present...skipping fltering step"
    fi
    
    ## Filtering the map file
    ##
    filteredGemStats=${filteredGem%.map.gz}.stats
    
    if [ $filteredGemStats -ot $filteredGem ];then
        step="GEM-STATS"
        startTime=$(date +%s)
        printHeader "Executing GEM stats step"
    
        ## Copy needed files to TMPDIR
        # copyToTmp "index"
    
        log "Producing stats for $filteredGem..." $step
        run "$gt_stats -i $filteredGem -t $threads -a -p 2> $filteredGemStats" "$ECHO"
        log "done\n" $step
        if [ -f $filteredGemStats ]; then
            log "Computing md5sum for stats file..." $step
            run "md5sum $filteredGemStats > $filteredGemStats.md5" "$ECHO"
            log "done\n"
        # else
        #    log "Error producing GEM stats" "ERROR" >&2
        #    exit -1
        fi
        endTime=$(date +%s)
        printHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "GEM stats file is present...skipping GEM stats step"
    fi
    
    ## Convert to bam and adding the XS field
    ##
    filteredBam=${sample}_filtered_cuff.bam
    indexFile="$BASEDIR/bp_rna_dashboard_20130809_temp.crg.txt"
    format="$BASEDIR/../tsv_format.json"
    
    if [ ! -e $filteredBam ];then
        step="CONVERT"
        startTime=$(date +%s)
        printHeader "Executing conversion step"
        
        ## Copy needed files to TMPDIR
        copyToTmp "index"
    
        log "Converting  $sample to bam..." $step
        readGroup=`$getHeaderMeta -i $indexFile -s $sample -f $format`
        if [[ $ECHO ]];then
            echo "$getHeaderMeta -i $indexFile -s $sample -f $format"
            echo "$readGroup"
        fi
        run "$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-33 -l --read-group $readGroup | sed 's/chrMT/chrM/g' | $addXS $readStrand | $samtools view -@ $threads -Sb - | $samtools sort -@ $threads -m 4G - $TMPDIR/${filteredBam%.bam}" "$ECHO"
        log "done\n" $step
        if [ -f $TMPDIR/$filteredBam ]; then
            log "Computing md5sum for filtered file..." $step
            run "md5sum $TMPDIR/$filteredBam > $TMPDIR/$filteredBam.md5" "$ECHO"
            run "cp $TMPDIR/$filteredBam.md5 ." "$ECHO"
            log "done\n"
            log "Copying filtered bam file to mapping dir..." $step
            run "cp $TMPDIR/$filteredBam ." "$ECHO"
            log "done\n"
        #else
        #    log "Error producing filtered bam file" "ERROR" >&2
        #    exit -1
        fi
        endTime=$(date +%s)
        printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Bam file is present...skipping conversion step"
    fi
else
    printHeader "Input file is $input...skipping mapping steps"
fi

#pipelineEnd=$(date +%s)
#    
#log "\n"
#printHeader "Blueprint pipeline for $sample completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
#exit 0

## Indexing the filtered bam file
##

[[ `basename $input` =~ bam ]] && filteredBam=`basename $input` && sample=`basename $input | sed 's/_filtered_cuff.bam//g'`

if [ $filteredBam.bai -ot $filteredBam ];then
    step="INDEX"
    startTime=$(date +%s)
    printHeader "Executing indexing step on the filtered bam file"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam"

    log "Indexing the filtered bam file\n" $step
    run "$samtools index $TMPDIR/$filteredBam" "$ECHO"
    if [ -f $TMPDIR/$filteredBam.bai ]; then
        log "Computing md5sum for filtered bam file index..." $step
        run "md5sum $TMPDIR/$filteredBam.bai > $TMPDIR/$filteredBam.bai.md5" "$ECHO"
        run "cp $TMPDIR/$filteredBam.bai.md5 ." "$ECHO"
        log "done\n"
        log "Copying bam index file to mapping dir..." $step
        run "cp $TMPDIR/$filteredBam.bai ." "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error producing bam index file" "ERROR" >&2
            exit -1
        fi
    fi
    endTime=$(date +%s)
    printHeader "Indexing step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Bam index file is present...skipping indexing step"
fi

## Producing stats for bam file
##
statsDir=$BASEDIR/$sample/stats
stats=true
if [ $stats ] && [ ! -d $statsDir ];then
    step="BAM-STATS"
    startTime=$(date +%s)
    printHeader "Executing bam stats step on the filtered bam file"

    if [ ! -d $statsDir ];then
        log "Creating stats folder..."
        run "mkdir $statsDir" "$ECHO"
        log "done\n"
    fi

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam"

    log "Producing mapping stats for the bam file\n" $step

    ## Create another bam with only uniquely mapping reads
    uniqBam=$TMPDIR/${filteredBam%.bam}_uniq.bam
    if [ ! -e $uniqBam ];then
        log "Making a bam file of unique mappings..." $step
        run "$BAMFLAG -in $filteredBam -out $uniqBam -m 3" "$ECHO"
        log "done\n"
    fi

    ## Create a bed12 from gtf
    genePred=$TMPDIR/${annName%.gtf}.GenePred
    bed12=$TMPDIR/${annName%.gtf}.bed12
    if [ ! -e $genePred ];then
        run "/users/rg/abreschi/bioprogs/Jim_Kent_source_tree/gtfToGenePred $annotation -allErrors $genePred 2> $genePred.err" "$ECHO"
        run "cat $genePred | ~abreschi/bioprogs/Jim_Kent_source_tree/genePredToBed12 > $bed12" "$ECHO"
    fi

    ## Check that the files were created successfully
    if [ ! -e $uniqBam ] || [ ! -e $genePred ];then
        if [[ ! $ECHO ]];then
            log "Error producing input file for stats" "ERROR" >&2
            exit -1
        fi
    fi

    ## Run the statistics
    #
    . /software/rg/el6.3/python2.7/bin/activate

    ## bam_stat
    bamStat=$statsDir/$sample.bam_stat
    if [ ! -e $bamStat ];then
        log "Computing general stats..." $step
        run "bam_stat.py -i $filteredBam 2> $bamStat" "$ECHO"
        log "done\n"
    fi


    ## Clipping_profile
    if [ `ls $statsDir/$sample.clipping_profile.* 2> /dev/null |wc -l` != 3 ];then
        log "Gettig clipping profile..." $step
        run "clipping_profile.py -i $filteredBam -o $statsDir/$sample" "$ECHO"
        log "done\n"
    fi

    # Gene body coverage
    # It outputs a file with counts and an R script. If the pdf is not produced,
    # the R script must be executed manually.
    # Stderr will contain a line like this:
    # Cannot generate pdf file from ERR180942_filtered_cuff.geneBodyCoverage_plot.r
    #
    if [ `ls $statsDir/$sample.geneBodyCoverage* 2> /dev/null |wc -l` != 3 ];then
        log "Gettig gene body coverage..." $step
        run "geneBody_coverage.py -r $bed12 -i $filteredBam -o $statsDir/$sample" "$ECHO"
        log "done\n"
    fi

    # Infer the configuration of the experiment, output is on stdout
    if [ ! -e $statsDir/$sample.inferred_exp ];then
        log "Inferring the configuration of the experiment..." $step
        run "infer_experiment.py -r $bed12 -i $filteredBam > $statsDir/$sample.inferred_exp" "$ECHO"
        log "done\n"
    fi

    # Insert size distribution
    if [ `ls $statsDir/$sample.inner_distance* 2> /dev/null |wc -l` != 4 ];then
        log "Getting the insert size distribution..." $step
        run "inner_distance.py -i $filteredBam -o $statsDir/$sample -r $bed12" "$ECHO"
        log "done\n"
    fi

    # Splice junctions
    if [ `ls $statsDir/$sample.[sj]* 2> /dev/null |wc -l` != 5 ];then
        log "Performing junction annotation..." $step
        run "junction_annotation.py -i $filteredBam -o $statsDir/$sample -r $bed12 &> $statsDir/$sample.junction_annotation.log" "$ECHO"
        log "done\n"
    fi

    run "deactivate" "$ECHO"

    endTime=$(date +%s)
    printHeader "BAM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Skipping BAM stats step"
fi

## Setting genome index files
#
genomeFai="$genome.fai"

## Producing bigWig files
##
doBigWig=0
if [[ $stranded == "1" ]];then
    eval "if [ ! -e $sample.plusRaw.bigwig ] || [ ! -e $sample.minusRaw.bigwig ];then doBigWig=1;fi"
else 
    eval "if [ ! -e $sample.bigwig ];then doBigWig=1;fi"
fi
if [[ $doBigWig == "1" ]];then
    step="BIGWIG"
    startTime=$(date +%s)
    printHeader "Executing BigWig step"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam"

    log "Producing bigWig files\n" $step

    if [[ $stranded == 1 ]];then

        ## Producing temporary bam with mate1 reversed
        revBam=$TMPDIR/${filteredBam%.bam}_rev1.bam
        log "Making temporary bam file with mate1 strand reversed..." $step
        run "$samtools view -h -@ $hthreads $TMPDIR/$filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $hthreads -Sb - > $revBam" "$ECHO"
        log "done\n"

        for strand in + -;
        do
            suffix="plusRaw"
            if [[ $strand == "-" ]];then
                suffix="minusRaw"
            fi
            bedGraph=$TMPDIR/$sample.$suffix.bedgraph
            bigWig=$TMPDIR/$sample.$suffix.bigwig
            log "Making bedGraph $strand strand\n" "$step"
            run "genomeCoverageBed -strand $strand -split -bg -ibam $revBam > $bedGraph" "$ECHO"
            log "Making bigWig $strand strand\n" $step
            run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"
        done

        if [ -f $TMPDIR/$sample.plusRaw.bigwig ] && [ -f $TMPDIR/$sample.minusRaw.bigwig ];then
            log "Computing md5sum for bigWig files..." $step
            run "md5sum $TMPDIR/$sample.plusRaw.bigwig > $TMPDIR/$sample.plusRaw.bigwig.md5" "$ECHO"
            run "md5sum $TMPDIR/$sample.minusRaw.bigwig > $TMPDIR/$sample.minusRaw.bigwig.md5" "$ECHO"
            run "cp $TMPDIR/$sample.plusRaw.bigwig.md5 ." "$ECHO"
            run "cp $TMPDIR/$sample.minusRaw.bigwig.md5 ." "$ECHO"
            log "done\n"
            log "Copying bigwig files to mapping dir..." $step
            run "cp $TMPDIR/$sample*.bigwig ." "$ECHO"
            log "done\n"
        else
            if [[ ! $ECHO ]];then
                log "Error producing bigWig files" "ERROR" >&2
                exit -1
            fi
        fi
    else
        bedGraph=$TMPDIR/$sample.bedgraph
        bigWig=$TMPDIR/$sample.bigwig
        log "Making bedGraph\n" "BEDGRAPH"
        run "genomeCoverageBed -split -bg -ibam $TMPDIR/${filteredBam} > $bedGraph" "$ECHO"
        log "Making bigWig\n" $step
        run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"

        if [ -f $bigWig ];then
            log "Computing md5sum for bigWig file..." $step
            run "md5sum $bigWig > $bigWig.md5" "$ECHO"
            run "cp $bigWig.md5 ." "$ECHO"
            log "done\n"
            log "Copying bigwig file to mapping dir..." $step
            run "cp $bigWig ." "$ECHO"
            log "done\n"
        else
            if [[ ! $ECHO ]];then
                log "Error producing bigWig files" "ERROR" >&2
                exit -1
            fi
        fi
    fi

    endTime=$(date +%s)
    printHeader "BigWig step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "BigWig files present...skipping bigwig step"
fi

## Setting genome index files
#
genomeFai="$genome.fai"

## Producing contig files
##
contigFile=${sample}_contigs.bed
if [ ! -e $contigFile ];then
    . /software/rg/el6.3/python2.7/bin/activate

    step="CONTIGS"
    startTime=$(date +%s)
    printHeader "Executing contigs step"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam"

    log "Producing contigs file\n" $step

    if [[ $stranded == 1 ]];then

        revBam=$TMPDIR/${filteredBam%.bam}_rev1.bam
        if [ ! -e $revBam ];then
            ## Producing temporary bam with mate1 reversed
            log "Making temporary bam file with mate1 strand reversed..." $step
            run "$samtools view -h -@ $hthreads $TMPDIR/$filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $hthreads -Sb - > $revBam" "$ECHO"
            log "done\n"
        fi

        uniqBam=${revBam%.bam}_uniq.bam
        if [ ! -e $uniqBam ];then
            log "Making a bam file of unique mappings..." $step
            run "$BAMFLAG -in $revBam -out $uniqBam -m 3" "$ECHO"
            log "done\n"
        fi


        for strand in + -;
        do
            suffix="plusRaw"
            if [[ $strand == "-" ]];then
                suffix="minusRaw"
            fi
            bedGraph=${uniqBam%.bam}.$suffix.bedgraph
            log "Making bedGraph $strand strand\n" "$step"
            run "genomeCoverageBed -strand $strand -split -bg -ibam $revBam > $bedGraph" "$ECHO"
        done

        log "Generationg the contigs file..." $step
        run "python $makecontig --chrFile $genomeFai --fileP ${uniqBam%.bam}.plusRaw.bedgraph --fileM ${uniqBam%.bam}.minusRaw.bedgraph | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\t\")} print s}' > $TMPDIR/$contigFile" "$ECHO"
        log "done\n"
    else
        
        uniqBam=$TMPDIR/${filteredBam%.bam}_uniq.bam
        if [ ! -e $uniqBam ];then
            log "Making a bam file of unique mappings..." $step
            run "$BAMFLAG -in $TMPDIR/$filteredBam -out $uniqBam -m 3" "$ECHO"
            log "done\n"
        fi
        
        log "Generationg the contigs file..." $step
        run "bamToBed -i $uniqBam | sort -k1,1 -nk2,2 | mergeBed > $TMPDIR/$contigFile" "$ECHO"
        log "done\n"
    fi

    if [ -f $TMPDIR/$contigFile ];then
        log "Computing md5sum for contigs file..." $step
        run "md5sum $TMPDIR/$contigFile > $TMPDIR/$contigFile.md5" "$ECHO"
        run "cp $TMPDIR/$contigFile.md5 ." "$ECHO"
        log "done\n"
        log "Copying contigs file to mapping dir..." $step
        run "cp $TMPDIR/$contigFile ." "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error producing contigs file" "ERROR" >&2
            exit -1
        fi
    fi

    run "deactivate" "$ECHO"

    endTime=$(date +%s)
    printHeader "Contigs step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Contigs file present...skipping contigs step"
fi

## Run FLux
##
step="FLUX"

quantDir="$BASEDIR/quantification"
export FLUX_MEM=16G

paramFile="$quantDir/param_file_Unstranded.par"
if [[ $stranded == 1 ]];then
    paramFile="$quantDir/param_file_Stranded.par"
fi

sample=${sample}
if [ ! -d $quantDir/$sample ]; then
    log "Creating sample folder in $quantDir..." $step
    run "mkdir -p $quantDir/$sample" "$ECHO"
    log "done\n"
fi

## Get profile file if required
#
if [[ $profile == 1 ]];then
    startTime=$(date +%s)
    printHeader "Executing Flux profiling step"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam,filtered-bai,annotation"

    if [ ! -e $TMPDIR/${annName%.gtf}_sorted.gtf ];then
        log "Checking if the annotation is sorted..." $step
        run "flux-capacitor -t sortGTF -c -i $TMPDIR/$annName -o $TMPDIR/${annName%.gtf}_sorted.gtf > $quantDir/$sample/${sample}_sort_annotation.log 2>&1" "$ECHO"
        log "done\n"
    fi

    annFile=$TMPDIR/${annName}
    if [ -e $TMPDIR/${annName%.gtf}_sorted.gtf ];then
        log "Using sorted annotation\n" $step
        annFile=$TMPDIR/${annName%.gtf}_sorted.gtf
    fi

    log "Preparing the parameter file..." $step
    run "cp $paramFile $TMPDIR" "$ECHO"
    paramFile=$TMPDIR/`basename $paramFile`
    run "echo \"PROFILE_FILE $TMPDIR/$sample.profile\" >> $paramFile" "$ECHO"
    #### TODO remove for new BAM files
    #run "echo \"USE_FLAGS false\" >> $paramFile" "$ECHO"
    log "done\n"

    if [ ! -e $quantDir/$sample/$sample.profile ];then
        log "Getting sistematic biases along transcripts\n" $step
        run "flux-capacitor --profile -p $paramFile -i $TMPDIR/$filteredBam -a $annFile > $quantDir/$sample/${sample}_flux_profile.log 2>&1" "$ECHO"
        if [ -e $TMPDIR/$sample.profile ];then
            log "Copying profile file to quantification dir..." $step
            run "cp $TMPDIR/$sample.profile $quantDir/$sample" "$ECHO"
        else
            if [[ ! $ECHO ]];then
                log "Error getting sistematic biases along transcripts" "ERROR" >&2
                exit -1
            fi
        fi
        log "done\n"
    else
        log "Copying profile file to $TMPDIR..." $step
        run "cp $quantDir/$sample/$sample.profile $TMPDIR" "$ECHO"
        log "done\n"
    fi
    endTime=$(date +%s)
    printHeader "Profiling step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
fi

## Run transcript quantification
#
if [ ! -e $quantDir/$sample/$sample.gtf ];then
    startTime=$(date +%s)
    printHeader "Executing Flux quantification step"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam,filtered-bai,annotation"

    if [ ! -e $TMPDIR/${annName%.gtf}_sorted.gtf ];then
        log "Checking if the annotation is sorted" $step
        run "flux-capacitor -t sortGTF -c -i $TMPDIR/$annName -o $TMPDIR/${annName%.gtf}_sorted.gtf > $quantDir/$sample/${sample}_sort_annotation.log 2>&1" "$ECHO"
        log "done\n"
    fi

    annFile=$TMPDIR/${annName}
    if [ -e $TMPDIR/${annName%.gtf}_sorted.gtf ];then
        log "Using sorted annotation\n" $step
        annFile=$TMPDIR/${annName%.gtf}_sorted.gtf
    fi

    if [[ $profile == 1 ]];then
        ## Copy profile file to TMPDIR if not there
        copyToTmp "flux-profile"
    fi

    if [ -e $TMPDIR/$sample.gtf ];then
        rm $TMPDIR/$sample.gtf
    fi

    log "Running Flux Capacitor\n" $step
    run "flux-capacitor -p $paramFile -i $TMPDIR/$filteredBam -a $TMPDIR/${annName} -o $TMPDIR/$sample.gtf > $quantDir/$sample/${sample}_flux_quantification.log 2>&1" "$ECHO"

    if [ -f $TMPDIR/$sample.gtf ]; then
        log "Computing md5sum for gtf file..." $step
        run "md5sum $TMPDIR/$sample.gtf > $TMPDIR/$sample.gtf.md5" "$ECHO"
        run "cp $TMPDIR/$sample.gtf.md5 $quantDir/$sample" "$ECHO"
        log "done\n"
        log "Copying gtf file to quantification dir..." $step
        run "cp $TMPDIR/$sample.gtf $quantDir/$sample" "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error running Flux Capacitor" "ERROR" >&2
            exit -1
        fi
    fi
    endTime=$(date +%s)
    printHeader "Quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Flux quantification file present...skipping Flux quantification step"
fi

txFile=$quantDir/$sample/${sample}_transcript.gtf
if [ ! -e $txFile ];then
    startTime=$(date +%s)
    printHeader "Getting transcript from quantifications"
    log "Generating transcripts file..."
    run "awk '\$3==\"transcript\"' $TMPDIR/$sample.gtf > $txFile" "$ECHO"
    log "done\n"
    log "Computing md5sum for transcripts file..." $step
    run "md5sum $txFile > $txFile.md5" "$ECHO"
    log "done\n"
    endTime=$(date +%s)
    printHeader "Transcripts written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Transcript file present...skipping transcript step"
fi

junctionFile=$quantDir/$sample/${sample}_junction.gtf
if [ ! -e $junctionFile ];then
    startTime=$(date +%s)
    printHeader "Getting junctions from quantifications"
    log "Generating junctions file..."
    run "awk '\$3==\"junction\"' $TMPDIR/$sample.gtf > $junctionFile" "$ECHO"
    log "done\n"
    log "Computing md5sum for junctions file..." $step
    run "md5sum $junctionFile > $junctionFile.md5" "$ECHO"
    log "done\n"
    endTime=$(date +%s)
    printHeader "Junctions written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Junctions file present...skipping transcript step"
fi

intronFile=$quantDir/$sample/${sample}_intron.gtf
if [ ! -e $intronFile ];then
    printHeader "Getting all-intronic regions from quantifications"
    log "Generating introns file..."
    run "awk '\$3==\"intron\"' $TMPDIR/$sample.gtf > $intronFile" "$ECHO"
    log "done\n"
    log "Computing md5sum for introns file..." $step
    run "md5sum $intronFile > $intronFile.md5" "$ECHO"
    log "done\n"
    endTime=$(date +%s)
    printHeader "Introns written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Introns file present...skipping transcript step"
fi

if [ ! -e $quantDir/$sample/${sample}_distinct_exon_with_rpkm.gff ];then
    step="EXON"
    startTime=$(date +%s)
    printHeader "Executing Exon quantification step"

    ## Copy needed files to TMPDIR
    copyToTmp "annotation,flux-gtf"

    log "Running Exon quantification\n" $step
    exd=$PWD
    run "cd $TMPDIR" "$ECHO"
    run "bash $trToEx $TMPDIR/$annName $TMPDIR/$sample.gtf" "$ECHO"
    run "cd $exd" "$ECHO"

    if [ -f $TMPDIR/${sample}_distinct_exon_with_rpkm.gff ]; then
        exonsFile="${sample}_distinct_exon_with_rpkm.gff"
        log "Computing md5sum for gff file..." $step
        run "md5sum $TMPDIR/$exonsFile > $TMPDIR/$exonsFile.md5" "$ECHO"
        cp $TMPDIR/$exonsFile.md5 $quantDir/$sample
        log "done\n"
        log "Copying gff file to quantification dir..." $step
        run "cp $TMPDIR/$exonsFile $quantDir/$sample" "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error running Exon quantification" "ERROR" >&2
            exit -1
        fi
    fi
    endTime=$(date +%s)
    printHeader "Exon quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

else
    printHeader "Exon quantification file present...skipping Exon quantification step"
fi

if [ ! -e $quantDir/$sample/${sample}_gene_with_rpkm.gff ];then
    step="GENE"
    startTime=$(date +%s)
    printHeader "Executing Gene quantification step"

    ## Copy needed files to TMPDIR
    copyToTmp "annotation,flux-gtf"

    log "Running Gene quantification\n" $step
    exd=$PWD
    run "cd $TMPDIR" "$ECHO"
    run "bash $trToGn $TMPDIR/$annName $TMPDIR/$sample.gtf" "$ECHO"
    run "cd $exd" "$ECHO"

    if [ -f $TMPDIR/${sample}_gene_with_rpkm.gff ]; then
        genesFile="${sample}_gene_with_rpkm.gff"
        log "Computing md5sum for gff file..." $step
        run "md5sum $TMPDIR/$genesFile > $TMPDIR/$genesFile.md5" "$ECHO"
        run "cp $TMPDIR/$genesFile.md5 $quantDir/$sample" "$ECHO"
        log "done\n"
        log "Copying gff file to quantification dir..." $step
        run "cp $TMPDIR/$genesFile $quantDir/$sample" "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error running Gene quantification" "ERROR" >&2
            exit -1
        fi
    fi
    endTime=$(date +%s)
    printHeader "Gene quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

else
    printHeader "Gene quantification file present...skipping Gene quantification step"
fi

pipelineEnd=$(date +%s)

log "\n"
printHeader "Blueprint pipeline for $sample completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
exit 0
