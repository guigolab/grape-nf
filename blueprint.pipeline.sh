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
    printf "  -t|--threads\t\tNumber of threads. Default \"1\".\n"
    printf "  -p|--paired-end\tSpecify whether the data is paired-end. Defalut: \"false\"\n"
    printf "  -c|--count-elements\tA comma separated list of elements to be counted by the Flux Capacitor.\n\t\t\tPossible values: INTRONS,SPLICE_JUNCTIONS. Defalut: \"none\"\n"
    printf "  -h|--help\t\tShow this message and exit.\n"
    printf "  --flux-mem\t\tSpecify the amount of ram the Flux Capacitor can use. Default: \"3G\".\n"
    printf "  --bam-stats\t\tRun the RSeQC stats on the bam file. Default \"false\".\n"
    printf "  --tmp-dir\t\tSpecify local temporary folder to copy files when running on distributed file systems. Default: \"-\".\n"
    printf "  --dry-run\t\tTest the pipeline. Writes the command to the standard output.\n"
    exit 0
}

function printHeader {
    local string=$1
    local message="*** $string ***"
    if [ ! $ECHO ]; then
        message="`date` $message"
    fi
    echo "$message"
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
    local dir=`dirname $path`
    local name=`basename $path`
    echo "`readlink -en $dir`/$name"
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
    paths=""
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        local name=`basename $i`
        local tmpfile="$tmpdir/$name"
        [[ $paths ]] && paths="$paths,"
        paths=${paths}${tmpfile}
        if [ ! -e $tmpfile ];then
            log "Copying $i to $TMPDIR..." $step
            run "cp $i $tmpdir" "$ECHO"
            log "done\n"
        else
            log "Skipping $i...temporary copy already exists\n" $step
        fi
    done
}

function finalizeStep {
    local file=$1
    local tmpdir=$2
    local outdir=$3
    paths=${outdir}/`basename $file`
    log "Computing md5sum for $file..." $step
    run "md5sum $file > $file.md5" "$ECHO"
    log "done\n"
    if [ -d $tmpdir ];then
        local tmpfile=$tmpdir/`basename $file`        
        log "Copying  temporary files..." $step
        run "cp $tmpfile $tmpfile.md5 $outdir" "$ECHO"
        log "done\n"
    fi
}

## Parsing arguments
#

# Execute getopt
ARGS=`getopt -o "i:g:a:m:n:s:t:l:q:r:c:hp" -l "input:,genome:,annotation:,mismatches:,hits:,read-strand:,threads:,loglevel:,quality:,max-read-length:,tmp-dir:,flux-mem:,count-elements:,bam-stats,dry-run,help,paired-end" \
      -n "$0" -- "$@"`

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
fluxMem="3G"
paired="false"

# general
loglevel="info"
threads="1"
tmpdir=${TMPDIR-"-"}
outdir=${SGE_O_WORKDIR-$PWD}

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
 
    -p|--paired-end)
      paired="true"
      shift ;;
 
    -c|--count-elements)
      if [ -n $2 ];
      then
        countElements="$2"
      fi
      shift 2;;

    --tmp-dir)
      if [ -n $2 ];
      then
        tmpdir=$(getAbsPath $2)
      fi
      shift 2;;
 
    --flux-mem)
      if [ -n $2 ];
      then
        fluxMem=$2
      fi
      shift 2;;

    --bam-stats)
      bamstats="true"
      shift ;;
    
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

BASEDIR=`dirname $outdir`
BINDIR="$BASEDIR/bin"

# activate python virtualenv
run ". $BASEDIR/bin/activate" "$ECHO"

export PATH=$BINDIR:$BINDIR/gemtools-1.6.2-i3/bin:$BINDIR/flux-capacitor-1.2.4/bin:$PATH

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
sample=${basename%.fastq*}
if [[ $paired == "true" ]];then
    sample=${basename%[_-\.]1*}
fi

genomeFai="$genome.fai"
gemIndex="${genome%.fa}.gem"
tindex="$annotation.junctions.gem"
tkeys="$annotation.junctions.keys"

annName=`basename $annotation`

hthreads=$((threads/2))
if [[ $hthreads == 0 ]];then
    hthreads=1
fi

## Binaries
#
gem2sam="gem-2-sam"
samtools="samtools"
trToGn="TrtoGn_RPKM.sh"
trToEx="TrtoEx_RPKM.sh"
bamToContigs="bamToContigs.sh"
gt_quality="gt.quality"
gt_filter="gt.filter"
gt_stats="gt.stats"
pigz="pigz"
bamflag="bamflag"
makecontig="contigsNew.py"
gtfToGenePred="gtfToGenePred"
genePredToBed12="genePredToBed12.awk"

## Print pipeline configuration

header="Pipeline configuration for $sample"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n" "Input file:" "$input"
printf "  %-34s %s\n" "Reference genome file:" "$genome"
printf "  %-34s %s\n" "Reference gene annotation file:" "$annotation"
printf "  %-34s %s\n" "Paired-end:" "$paired"
printf "  %-34s %s\n" "Max number of allowed mismatches:" "$mism"
printf "  %-34s %s\n" "Max number of hits:" "$hits"
printf "  %-34s %s\n" "Quality offset:" "$qualityOffset"
printf "  %-34s %s\n" "Max read length:" "$maxReadLength"
printf "  %-34s %s\n" "Strandedness:" "$readStrand"
printf "  %-34s %s\n" "Number of threads:" "$threads"
printf "  %-34s %s\n" "Flux Capacitor memory:" "$fluxMem"
printf "  %-34s %s\n" "Temporary folder:" "$tmpdir"
printf "  %-34s %s\n" "Loglevel:" "$loglevel"
printf "\n\n"

## START
#

printHeader "Starting Blueprint pipeline for $sample"
pipelineStart=$(date +%s)

## Mapping
#
if [[ `basename $input` =~ fastq ]];then 
    gem="$outdir/$sample.map.gz"
    if [ ! -e $gem ];then
        step="MAP"
        startTime=$(date +%s)
        printHeader "Executing mapping step"                
            
        if [ -d $tmpdir ]; then
            ## Copy needed files to TMPDIR
            copyToTmp "$gemIndex,$annotation,$tindex,$tkeys"
            IFS=',' read index annotation tindex tkeys <<< "$paths"
            gem="$tmpdir/$sample.map.gz"
        fi
    
        log "Running gemtools rna pipeline on ${sample}" $step
        command="gemtools --loglevel $loglevel rna-pipeline -f $input -q $qualityOffset -i $gemIndex -a $annotation -o `dirname $gem` -t $threads --no-stats --no-bam"
        if [[ $paired == "false" ]]; then
            command="$command --single-end"
        fi
        run "$command" "$ECHO"
   
        set -e && finalizeStep $gem $tmpdir $outdir
        IFS=',' read gem <<< "$paths"
        
        endTime=$(date +%s)
        printHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Map file already present...skipping mapping step"
    fi
 
    ## Filtering the map file
    ##
    filteredGem=${gem%.map.gz}_m${mism}_n${hits}.map.gz    
    
    if [ ! -e $filteredGem ];then
        step="FILTER"
        startTime=$(date +%s)
        printHeader "Executing filtering step"
        
        log "Filtering map file..." $step
        run "$gt_quality -i $gem -t $threads | $gt_filter --max-levenshtein-error $mism -t $threads | $gt_filter --max-matches $hits -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
        log "done\n" $step

        set -e && finalizeStep $filteredGem "-" $outdir

        endTime=$(date +%s)
        printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Filtered map file is present...skipping fltering step"
    fi
 
    ## Getting stats for the filtered map file
    ##
    filteredGemStats=${filteredGem%.map.gz}.stats
    
    if [ $filteredGemStats -ot $filteredGem ]; then
        step="GEM-STATS"
        startTime=$(date +%s)
        printHeader "Executing GEM stats step"
    
        log "Producing stats for $filteredGem..." $step
        run "$gt_stats -i $filteredGem -t $threads -a -p 2> $filteredGemStats" "$ECHO"
        log "done\n" $step

        set -e && finalizeStep $filteredGemStats "-" $outdir

        endTime=$(date +%s)
        printHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "GEM stats file is present...skipping GEM stats step"
    fi

    ## Convert to bam 
    ##
    filteredBam=${filteredGem%.map.gz}.bam
    
    if [ ! -e $filteredBam ]; then
        step="CONVERT"
        startTime=$(date +%s)
        printHeader "Executing conversion step"

        if [ -d $tmpdir ]; then
            ## Copy needed files to TMPDIR
            copyToTmp "$gemIndex"
            IFS=',' read index <<< "$paths"
            filteredBam="$tmpdir/`basename $filteredBam`"
        fi
       
        log "Converting  $sample to bam..." $step
        command="$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $gemIndex -q offset-$qualityOffset -l"
        if [[ $paired == "true" ]]; then
            command="$command --expect-paired-end-reads"
        else
            command="$command --expect-single-end-reads | awk 'BEGIN{OFS=FS=\"\t\"}\$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and(\$2,a[i])>0){\$2=xor(\$2,a[i])}}}{print}'"
        fi

        if [[ $readGroup ]];
        then
            command="$command --read-group $readGroup"
        fi
        run "$command | $samtools view -@ $threads -Sb - | $samtools sort -@ $threads -m 4G - ${filteredBam%.bam}" "$ECHO"
        log "done\n" $step
        
        set -e && finalizeStep $filteredBam $tmpdir $outdir       
        IFS=',' read filteredBam <<< "$paths"

        endTime=$(date +%s)
        printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Bam file is present...skipping conversion step"
    fi
else
    printHeader "Input file is $input...skipping mapping steps"
fi

## Indexing the filtered bam file
##
filteredBai="$filteredBam.bai"
if [ $filteredBai -ot $filteredBam ];then
    step="INDEX"
    startTime=$(date +%s)
    printHeader "Executing indexing step on the filtered bam file"

    if [ -d $tmpdir ]; then
        ## Copy needed files to TMPDIR
        copyToTmp "$filteredBam"
        IFS=',' read filteredBam <<< "$paths"
        filteredBai="$tmpdir/`basename $filteredBai`"
    fi

    log "Indexing the filtered bam file\n" $step
    run "$samtools index $filteredBam" "$ECHO"
    
    set -e && finalizeStep $filteredBai $tmpdir $outdir        
    IFS=',' read filteredBai <<< "$paths"

    endTime=$(date +%s)
    printHeader "Indexing step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Bam index file is present...skipping indexing step"
fi

## Producing stats for bam file
##
statsDir=$outdir/stats
if [ $bamstats ]; then
    step="BAM-STATS"
    startTime=$(date +%s)
    printHeader "Executing bam stats step on the filtered bam file"

    if [ ! -d $statsDir ];then
        log "Creating stats folder..."
        run "mkdir $statsDir" "$ECHO"
        log "done\n"
    fi

    if [ -d $tmpdir ]; then
        ## Copy needed files to TMPDIR
        copyToTmp "$filteredBam,$annotation"
        IFS=',' read filteredBam annotation <<< "$paths"
    fi

    log "Producing mapping stats for the bam file\n" $step

    ## Create another bam with only uniquely mapping reads
    uniqBam=${filteredBam%.bam}_uniq.bam
    if [ ! -e $uniqBam ];then
        log "Making a bam file of unique mappings..." $step
        set +e && run "$bamflag -in $filteredBam -out $uniqBam -m 3" "$ECHO"
        log "done\n"
    fi

    ## Create a bed12 from gtf
    genePred=${annotation%.gtf}.GenePred
    bed12=${annotation%.gtf}.bed12
    if [ ! -e $genePred ];then
        set +e && run "$gtfToGenePred $annotation -allErrors $genePred 2> $genePred.err" "$ECHO"
        set +e && run "cat $genePred | $genePredToBed12 > $bed12" "$ECHO"
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

    # Splice junctions
    if [ `ls $statsDir/$sample.[sj]* 2> /dev/null |wc -l` != 5 ];then
        log "Performing junction annotation..." $step
        run "junction_annotation.py -i $filteredBam -o $statsDir/$sample -r $bed12 &> $statsDir/$sample.junction_annotation.log" "$ECHO"
        log "done\n"
    fi

    if [[ $paired == "true" ]]; then
        # Insert size distribution
        if [ `ls $statsDir/$sample.inner_distance* 2> /dev/null |wc -l` != 4 ];then
            log "Getting the insert size distribution..." $step
            run "inner_distance.py -i $filteredBam -o $statsDir/$sample -r $bed12" "$ECHO"
            log "done\n"
        fi
    fi

    endTime=$(date +%s)
    printHeader "BAM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Skipping BAM stats step"
fi


## Producing bigWig files
##
doBigWig=0
if [[ $readStrand != "NONE" ]];then
    eval "if [ ! -e $sample.plusRaw.bigwig ] || [ ! -e $sample.minusRaw.bigwig ];then doBigWig=1;fi"
else 
    eval "if [ ! -e $sample.bigwig ];then doBigWig=1;fi"
fi
if [[ $doBigWig == "1" ]];then
    step="BIGWIG"
    startTime=$(date +%s)
    printHeader "Executing BigWig step"

    ## Copy needed files to TMPDIR
    if [ -d $tmpdir ]; then
        copyToTmp "$filteredBam"
        IFS=',' read filteredBam <<< "$paths"
    fi

    log "Producing bigWig files\n" $step

    if [[ $readStrand != "NONE" ]];then

        ## Producing temporary bam with mate1 reversed
        revBam="${filteredBam%.bam}_1rev.bam"
        log "Making temporary bam file with mate1 strand reversed..." $step
        run "$samtools view -h -@ $hthreads $filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $hthreads -Sb - > $revBam" "$ECHO"
        log "done\n"

        for strand in + -;
        do
            suffix="plusRaw"
            if [[ $strand == "-" ]];then
                suffix="minusRaw"
            fi

            bedGraph=$outdir/$sample.$suffix.bedgraph
            bigWig=$outdir/$sample.$suffix.bigwig
            if [ -d $tmpdir ]; then                
                bedGraph=$tmpdir/$sample.$suffix.bedgraph
                bigWig=$tmpdir/$sample.$suffix.bigwig
            fi
            log "Making bedGraph $strand strand\n" "$step"
            run "genomeCoverageBed -strand $strand -split -bg -ibam $revBam > $bedGraph" "$ECHO"
            log "Making bigWig $strand strand\n" $step
            run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"

            set -e && finalizeStep $bigWig $tmpdir $outdir
            IFS=',' read bigWig <<< "$paths"
        done
    else
        bedGraph=$outdir/$sample.bedgraph
        bigWig=$outdir/$sample.bigwig
        if [ -d $tmpdir ]; then                
            bedGraph=$tmpdir/$sample.bedgraph
            bigWig=$tmpdir/$sample.bigwig
        fi

        log "Making bedGraph\n" "BEDGRAPH"
        run "genomeCoverageBed -split -bg -ibam $filteredBam > $bedGraph" "$ECHO"
        log "Making bigWig\n" $step
        run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"

        set -e && finalizeStep $bigWig $tmpdir $outdir
        IFS=',' read bigWig <<< "$paths"
    fi

    endTime=$(date +%s)
    printHeader "BigWig step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "BigWig files present...skipping bigwig step"
fi

## Producing contig files
##
contigFile=$outdir/${sample}_contigs.bed
if [ ! -e $contigFile ];then
    step="CONTIGS"
    startTime=$(date +%s)
    printHeader "Executing contigs step"

    if [ -d $tmpdir ]; then
        ## Copy needed files to TMPDIR
        copyToTmp "$filteredBam"
        IFS=',' read filteredBam <<< "$paths"
        contigFile=$tmpdir/${sample}_contigs.bed
    fi

    log "Producing contigs file\n" $step

    if [[ $readStrand != "NONE" ]];then

        revBam=${filteredBam%.bam}_1rev.bam
        if [ ! -e $revBam ];then
            ## Producing temporary bam with mate1 reversed
            log "Making temporary bam file with mate1 strand reversed..." $step
            run "$samtools view -h -@ $threads $filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $threads -Sb - > $revBam" "$ECHO"
            log "done\n"
        fi

        uniqBam=${revBam%.bam}_uniq.bam
        if [ ! -e $uniqBam ];then
            log "Making a bam file of unique mappings..." $step
            set +e && run "$bamflag -in $revBam -out $uniqBam -m 3" "$ECHO"
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
        run "python $makecontig --chrFile $genomeFai --fileP ${uniqBam%.bam}.plusRaw.bedgraph --fileM ${uniqBam%.bam}.minusRaw.bedgraph | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\t\")} print s}' > $contigFile" "$ECHO"
        log "done\n"
    else
        
        uniqBam=${filteredBam%.bam}_uniq.bam
        if [ ! -e $uniqBam ];then
            log "Making a bam file of unique mappings..." $step
            set +e && run "$bamflag -in $filteredBam -out $uniqBam -m 3" "$ECHO"
            log "done\n"
        fi
        
        log "Generationg the contigs file..." $step
        run "bamToBed -i $uniqBam | sort -k1,1 -nk2,2 | mergeBed > $contigFile" "$ECHO"
        log "done\n"
    fi

    set -e && finalizeStep $contigFile $tmpdir $outdir
    IFS=',' read contigFile <<< "$paths"

    endTime=$(date +%s)
    printHeader "Contigs step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Contigs file present...skipping contigs step"
fi

## Run FLux
##
step="FLUX"

quantDir="$BASEDIR/quantification/$sample"
export FLUX_MEM=$fluxMem

if [ ! -d $quantDir ]; then
    log "Creating sample folder in $quantDir..." $step
    run "mkdir -p $quantDir" "$ECHO"
    log "done\n"
fi

paramFile="$quantDir/$sample.par"
proFile="$quantDir/$sample.profile"

# prepare parameter file
#
# READ_STRAND MATE2_SENSE
# ANNOTATION_MAPPING PAIRED_STRANDED
# COUNT_ELEMENTS [SPLICE_JUNCTIONS, INTRONS]

if [ ! -e $paramFile ]; then
    run "echo \"# Flux Capacitor parameter file for $sample\" >> $paramFile" "$ECHO"
    annotationMapping="AUTO"
    if [[ $readStrand != "NONE" ]]; then
        run "echo \"READ_STRAND $readStrand\" >> $paramFile" "$ECHO"
        annotationMapping="STRANDED"
        if [[ $paired == "true" ]];then 
            annotationMapping="PAIRED_${annotationMapping}"
        else
            annotationMapping="SINGLE_${annotationMapping}"
        fi
    fi
    
    run "echo \"ANNOTATION_MAPPING $annotationMapping\" >> $paramFile" "$ECHO"    
    run "echo \"COUNT_ELEMENTS [$countElements]\" >> $paramFile" "$ECHO"
fi

## Show Flux parameter file
#
run "echo \"\"" "$ECHO"
run "cat $paramFile" "$ECHO"
run "echo \"\"" "$ECHO"

## Run transcript quantification
#
fluxGtf="$quantDir/$sample.gtf"
if [ ! -e $fluxGtf ];then
    startTime=$(date +%s)
    printHeader "Executing Flux quantification step"

    if [ -d $tmpdir ];then
        ## Copy needed files to TMPDIR
        copyToTmp "$filteredBam,$filteredBai,$annotation"
        IFS=',' read filteredBam filteredBai annotation <<< "$paths"
        fluxGtf="$tmpdir/$sample.gtf"
    fi

    if [ ! -e $TMPDIR/${annName%.gtf}_sorted.gtf ];then
        sortLog="$quantDir/${sample}_sort_annotation.log"
        log "Checking if the annotation is sorted" $step
        set -e && run "flux-capacitor -t sortGTF -c -i $annotation -o ${annotation%.gtf}_sorted.gtf > $sortLog 2>&1" "$ECHO"
        log "done\n"
    fi

    anno=$annotation
    if [ -e {$anno%.gtf}_sorted.gtf ];then
        anno=${anno%.gtf}_sorted.gtf
        log "Using sorted annotation: $anno\n" $step
    fi

    if [[ ! -e $proFile ]];then
        profileLog="$quantDir/${sample}_flux_profile.log"
        log "Getting sistematic biases along transcripts\n" $step
        run "flux-capacitor --profile -p $paramFile -i $filteredBam -a $anno --profile-file $proFile > $profileLog 2>&1" "$ECHO"
    fi

    if [ -e $fluxGtf ];then
        rm $fluxGtf
    fi

    log "Running Flux Capacitor\n" $step
    quantLog=$quantDir/${sample}_flux_quantification.log
    run "flux-capacitor -p $paramFile -i $filteredBam -a $anno -o $fluxGtf > $quantLog 2>&1" "$ECHO"

    set -e && finalizeStep $fluxGtf $tmpdir $quantDir
    IFS=',' read fluxGtf <<< "$paths"

    endTime=$(date +%s)
    printHeader "Quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Flux quantification file present...skipping Flux quantification step"
fi

txFile=$quantDir/${sample}_transcript.gtf
if [ ! -e $txFile ];then
    startTime=$(date +%s)
    printHeader "Getting transcript from quantifications"
    log "Generating transcripts file..."
    run "awk '\$3==\"transcript\"' $fluxGtf > $txFile" "$ECHO"
    log "done\n"
    log "Computing md5sum for transcripts file..." $step
    run "md5sum $txFile > $txFile.md5" "$ECHO"
    log "done\n"
    endTime=$(date +%s)
    printHeader "Transcripts written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Transcript file present...skipping transcript step"
fi

if [[ $countElements =~ "SPLICE_JUNCTIONS" ]]; then
    junctionFile=$quantDir/${sample}_junction.gtf
    if [ ! -e $junctionFile ];then
        startTime=$(date +%s)
        printHeader "Getting junctions from quantifications"
        log "Generating junctions file..."
        run "awk '\$3==\"junction\"' $fluxGtf > $junctionFile" "$ECHO"
        log "done\n"
        log "Computing md5sum for junctions file..." $step
        run "md5sum $junctionFile > $junctionFile.md5" "$ECHO"
        log "done\n"
        endTime=$(date +%s)
        printHeader "Junctions written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Junctions file present...skipping transcript step"
    fi
fi

if [[ $countElements =~ "INTRONS" ]]; then
    intronFile=$quantDir/${sample}_intron.gtf
    if [ ! -e $intronFile ];then
        printHeader "Getting all-intronic regions from quantifications"
        log "Generating introns file..."
        run "awk '\$3==\"intron\"' $fluxGtf > $intronFile" "$ECHO"
        log "done\n"
        log "Computing md5sum for introns file..." $step
        run "md5sum $intronFile > $intronFile.md5" "$ECHO"
        log "done\n"
        endTime=$(date +%s)
        printHeader "Introns written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else
        printHeader "Introns file present...skipping transcript step"
    fi
fi

exonFile=$quantDir/${sample}_distinct_exon_with_rpkm.gff
if [ ! -e $exonFile ];then
    step="EXON"
    startTime=$(date +%s)
    printHeader "Executing Exon quantification step"

    
    if [ -d $tmpdir ]; then
        ## Copy needed files to TMPDIR
        copyToTmp "$annotation,$fluxGtf"
        IFS=',' read annotation fluxGtf <<< "$paths"
        exonFile=$tmpdir/${sample}_distinct_exon_with_rpkm.gff
    fi

    log "Running Exon quantification\n" $step
    run "$trToEx -a $annotation -i $fluxGtf -o `dirname $exonFile`" "$ECHO"

    set -e && finalizeStep $exonFile $tmpdir $quantDir
    IFS=',' read exonFile <<< "$paths"
    
    endTime=$(date +%s)
    printHeader "Exon quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

else
    printHeader "Exon quantification file present...skipping Exon quantification step"
fi

geneFile=$quantDir/${sample}_gene_with_rpkm.gff
if [ ! -e $geneFile ];then
    step="GENE"
    startTime=$(date +%s)
    printHeader "Executing Gene quantification step"

    if [ -d $tmpdir ]; then
        ## Copy needed files to TMPDIR
        copyToTmp "$annotation,$fluxGtf"
        IFS=',' read annotation fluxGtf <<< "$paths"
        geneFile=$tmpdir/${sample}_gene_with_rpkm.gff
    fi

    log "Running Gene quantification\n" $step
    run "$trToGn -a $annotation -i $fluxGtf -o `dirname $geneFile`" "$ECHO"

    set -e && finalizeStep $geneFile $tmpdir $quantDir
    IFS=',' read geneFile <<< "$paths"
    
    endTime=$(date +%s)
    printHeader "Gene quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

else
    printHeader "Gene quantification file present...skipping Gene quantification step"
fi

# deactivate python virtualenv
run "deactivate" "$ECHO"

pipelineEnd=$(date +%s)

log "\n"
printHeader "Blueprint pipeline for $sample completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "

exit 0
