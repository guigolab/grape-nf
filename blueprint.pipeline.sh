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
#$ -q rg-el6
#$ -l virtual_free=64G
#
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err
#
function usage {
    echo "Usage: $0 -i <fastq_file> -s <sex> [OPTION]..."
    echo "Execute the Blueprint pipeline on one sample."
    echo ""
    printf "\t-i\tinput file\n"
    printf "\t-g\tspecify the sex of the sample [M | F]\n"
    echo ""
    echo "Options:"
    printf "\t-s\tflag to specify whether the sample has strandness information for the reads. Deafault \"false\"\n"
    printf "\t-d\tdirectionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default \"NONE\".\n"
    printf "\t-l\tLog level (error, warn, info, debug). Default \"info\".\n"
    printf "\t-p\twrite Flux Capacitor profile file. Default \"false\".\n"
    printf "\t-t\tTest the pipeline. Writes the command to the standard output.\n"
    exit 0
}

function printHeader {
    string=$1
    echo "`date` *** $string ***"
}

function log {
    string=$1
    label=$2
    if [[ ! $ECHO ]];then
        if [[ $label != "" ]];then
            printf "[$label] $string"
        else
            printf "$string"
        fi
    fi
}

function run {
    command=($1)
    if [[ $2 ]];then
         ${2}${command[@]}
    else
        eval ${command[@]}
    fi

}

function copyToTmp {
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        case $i in
            "annotation")
                if [ ! -e $TMPDIR/$annName ];then
                    log "Copying annotation file to $TMPDIR..." $step
                    run "cp $annotation $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "index")
                if [ ! -e $TMPDIR/`basename $index` ];then
                    log "Copying genome index file to $TMPDIR..." $step
                    run "cp $index $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "t-index")
                if [ ! -e $TMPDIR/$annName.gem ];then
                    log "Copying annotated transcriptome index file to $TMPDIR..." $step
                    run "cp $annotation.gem $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "keys")
                if [ ! -e $TMPDIR/$annName.junctions.keys ];then
                    log "Copying annotated transcriptome keys file to $TMPDIR..." $step
                    run "cp $annotation.junctions.keys $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "fastq")
                if [ ! -e $TMPDIR/`basename $input` ];then
                    log "Copying fastq files to $TMPDIR..." $step
                    run "cp $input $TMPDIR" "$ECHO"
                    run "cp ${input/_1.fastq/_2.fastq} $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "map.gz")
                if [ ! -e $TMPDIR/$sample.map.gz ]; then
                    log "Copying map file to $TMPDIR..." $step
                    run "cp $sample.map.gz $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "bam")
                if [ ! -e $TMPDIR/$sample.bam ]; then
                    log "Copying bam file to $TMPDIR..." $step
                    run "cp $sample.bam $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "filtered-bam")
                if [ ! -e $TMPDIR/$filteredBam ]; then
                    log "Copying filtered bam file to $TMPDIR..." $step
                    run "cp $filteredBam $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "filtered-bai")
                if [ ! -e $TMPDIR/$filteredBam.bai ];then
                    log "Copying index for filtered bam file to $TMPDIR..." $step
                    run "cp $filteredBam.bai $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "flux-profile")
                if [ ! -e $TMPDIR/$sample.profile ];then
                    log "Copying profile file to $TMPDIR..." $step
                    run "cp $quantDir/$sample/$sample.profile $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "flux-gtf")
               if [ ! -e $TMPDIR/$sample.gtf ];then
                   log "Copying Flux file to $TMPDIR..." $step
                   run "cp $quantDir/$sample/$sample.gtf $TMPDIR" "$ECHO"
                   log "done\n"
               fi
            esac
    done
}

## Parsing arguments
#

while getopts ":i:g:std:l:ph" opt; do
  case $opt in
    i)
      input=$OPTARG
      ;;
    g)
      sex=$OPTARG
      ;;
    s)
      stranded=1
      ;;
    t)
      ECHO="echo "
      ;;
    d)
      readStrand=$OPTARG
      ;;
    l)
      loglevel=$OPTARG
      ;;
    p)
      profile=1
      ;;
    h)
      usage
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

## Setting up the environment
#
BASEDIR=`dirname ${SGE_O_WORKDIR-$PWD}`
export PATH=$BASEDIR/gemtools-1.6-i3/bin:/software/rg/el6.3/flux-capacitor-1.2.4-SNAPSHOT/bin:$HOME/bin:$PATH

## Setting variables and input files
##
if [[ $input == "" ]];then
    log "Please specify the input file\n" "ERROR" >&2
    exit -1
fi

if [[ $sex == "" ]] || [[ $sex != [MFU] ]];then
    log "Please specify the sex\n" "ERROR" >&2
    exit -1
fi

if [[ $readStrand == "" ]];then
    readStrand="NONE"
fi

if [[ $loglevel == "" ]];then
    loglevel="info"
fi

basename=$(basename $input)
sample=${basename%_1*}
threads=${NSLOTS-1}

index="$BASEDIR/Homo_sapiens.GRCh37.chromosomes.chr.gem"
annotation="$BASEDIR/gencode.v15.annotation.gtf"
if [[ $sex == "F" ]];then
    index="$BASEDIR/Homo_sapiens.GRCh37.chromosomes.female.chr.gem"
    annotation="$BASEDIR/gencode.v15.annotation.female.gtf"
fi

annName=`basename $annotation`

gem2sam="/users/rg/epalumbo/projects/GTEx/gem.evaluation/mapping/convert/gem-2-sam"
#samtools="/users/rg/epalumbo/projects/GTEx/gem.evaluation/mapping/convert/samtools -n $threads"
samtools="/users/rg/epalumbo/git/samtools/samtools"
addXS="/users/rg/abreschi/Documents/utils/sam2cufflinks.sh"
trToGn="/users/rg/sdjebali/bin/TrtoGn_RPKM.sh"
trToEx="/users/rg/sdjebali/bin/TrtoEx_RPKM.sh"
bamToContigs="~sdjebali/bin/TrtoEx_RPKM.sh"
filter="$BASEDIR/bin/filterMapq78.sh"

## START
##
printHeader "Starting Blueprint pipeline for $sample"
pipelineStart=$(date +%s)

## Mapping
##
if [ ! -e $sample.map.gz ];then
    step="MAP"
    startTime=$(date +%s)
    printHeader "Executing mapping step"

    ## Activate the python virtualenv
    run ". $BASEDIR/venv/bin/activate" "$ECHO"

    ## Copy needed files to TMPDIR
    copyToTmp "index,annotation,t-index,keys"

    log "Running gemtools rna pipeline on ${sample}" $step
    run "gemtools --loglevel $loglevel rna-pipeline -f $input -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -t $threads --no-bam" "$ECHO"
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

## Converting to bam
##

hthreads=$((threads/2))
if [[ $hthreads == 0 ]];then
    hthreads=1
fi


if [ ! -e $sample.bam ];then
    step="CONVERT"
    startTime=$(date +%s)
    printHeader "Executing conversion step"

    ## Copy needed files to TMPDIR
    copyToTmp "index"

    log "Converting ${sample} to bam\n" $step

    run "pigz -p $hthreads -dc $sample.map.gz | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-33 | $samtools view -@ $hthreads -Sb - | $samtools sort -@ $hthreads -m `echo $((4<<30))` - $sample" "$ECHO"
    #pigz -p $threads -dc $TMPDIR/$sample.map.gz | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-33 | $samtools view -Sb - | $samtools sort -m $((8<<30)) - $TMPDIR/$sample
    if [ -f $TMPDIR/${sample}.bam ]; then
        log "Computing md5sum for bam file..." $step
        run "md5sum $TMPDIR/$sample.bam > $TMPDIR/$sample.bam.md5" "$ECHO"
        run "cp $TMPDIR/$sample.bam.md5 ." "$ECHO"
        log "done\n"

        log "Copying bam file to mapping dir..." $step
        run "cp $TMPDIR/${sample}.bam ." "$ECHO"
        log "done\n"
    #else
    #    log "Error producing bam file" "ERROR" >&2
    #    exit -1
    fi
    endTime=$(date +%s)
    printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Bam file already present...skipping conversion step"
fi

## Filtering the bam and adding the XS field
##
filteredBam=${sample}_filtered_cuff.bam

if [ ! -e $filteredBam ];then
    step="FILTER"
    startTime=$(date +%s)
    printHeader "Executing filtering step"

    log "Filtering bam file and adding XS field\n" $step
    run "$samtools view -@ $hthreads -h $sample.bam |  awk -F \"\t\" '\$1~/@/ || and(\$2,4)>0 || (\$5>=78 && (\$5<=90 || \$5>=119))' | sed 's/chrMT/chrM/g' | $addXS $readStrand | $samtools view -@ $hthreads -Sb - > $filteredBam" "$ECHO"
    #$samtools view -h $TMPDIR/$sample.bam | awk -F"\t" '$1~/@/ || and($2,4)>0 || ($5>=78 && ($5<=90 || $5>=119))' | $addXS $readStrand | $samtools view -Sb - > $TMPDIR/$filteredBam
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
    printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Filtered bam file is present...skipping fltering step"
fi

## Indexing the filtered bam file
##

if [ ! -e $filteredBam.bai ];then
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
        log "Copying filtered bam file index to mapping dir..." $step
        run "cp $TMPDIR/$filteredBam.bai ." "$ECHO"
        log "done\n"
    else
        if [[ ! $ECHO ]];then
            log "Error producing filtered bam file index" "ERROR" >&2
            exit -1
        fi
    fi
    endTime=$(date +%s)
    printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Index for fitlered bam file is present...skipping indexing step"
fi


## Producing bigWig files
##

if [ ! -e $sample.plusRaw.bigwig ] || [ ! -e $sample.minusRaw.bigwig ];then
    step="BIGWIG"
    startTime=$(date +%s)
    printHeader "Executing BigWig step"

    ## Copy needed files to TMPDIR
    copyToTmp "filtered-bam"

    log "Producing bigWig files\n" $step

    indexDir="/users/rg/projects/references/Genome"
    genomeFai="$indexDir/Homo_sapiens.GRCh37.chromosomes.chr.M.fa.fai"
    if [[ $sex == "F" ]];then
        genomeFai="$indexDir/Homo_sapiens.GRCh37.chromosomes.female.chr.M.fa.fai"
    fi

    for strand in + -;
    do
        suffix="plusRaw"
        if [[ $strand == "-" ]];then
            suffix="minusRaw"
        fi
        bedGraph=$TMPDIR/$sample.$suffix.bedgraph
        bigWig=$TMPDIR/$sample.$suffix.bigwig
        log "Making bedGraph $strand strand\n" "BEDGRAPH"
        run "genomeCoverageBed -strand $strand -split -bg -ibam $TMPDIR/$filteredBam > $bedGraph" "$ECHO"
        log "Making bigWig $strand strand\n" $step
        run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"
    done

    if [ -f $TMPDIR/$sample.plusRaw.bigwig ] && [ -f $TMPDIR/$sample.minusRaw.bigwig ];then
        log "Computin md5sum for bigWig files..." $step
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
    endTime=$(date +%s)
    printHeader "BigWig step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "BigWig files present...skipping bigwig step"
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
