#!/usr/bin/env nextflow

/*
 * Copyright (c) 2013, Centre for Genomic Regulation (CRG)
 * Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
 *
 * This file is part of the Blueprint RNAseq pipeline.
 *
 * The Blueprint RNAseq pipeline is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

env = System.getenv()

params.input = 'test/*.fastq.gz'                               //The input file
params.genome = 'tutorial/data/genome_1Mbp.fa'                              //The reference genome file
params.annotation = 'tutorial/data/annotation.gtf'                          //The reference gene annotation file
params.mismatches = 4                           //Max number of mismatches. Default 4	
params.hits = 10                                //Max number of hits. Default 10
params.quality_offset = 33                      //The quality offset of the fastq files. Default: 33
params.max_read_length = 150                    //The maximum read length (used to compute the transcriptomes). Default: 150
params.read_strand = 'NONE'                     //The directionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default NONE
params.loglevel = 'WARN'                        //Log level (error, warn, info, debug). Default info
params.cpus = 1                                 //Number of threads. Default 1
params.paired_end = false                       //Specify whether the data is paired-end. Default: false 
params.count_elements = []                      //A comma separated list of elements to be counted by the Flux Capacitor. Possible values: INTRONS,SPLICE_JUNCTIONS. Defalut: none	
params.read_group = ''                          //A comma separated list of tags for the @RG field of the BAM file. Check the SAM specification for details. Default: none
params.bam_stats = false                        //Run the RSeQC stats on the bam file. Default false
params.flux_mem	= '3G'                          //The amount of ram the Flux Capacitor can use. Default: 3G
params.tmp_dir = (env.TMPDIR != null ? true : false)                             //The local temporary folder to copy files when running on shared file systems. Default: TMPDIR
params.outdir = "$PWD"                          //The general output folder
params.steps = 'mapping,bigwig,contig,flux'     //The steps to be executed

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

// Setting up environment
baseDir = file(config.baseDir)
binDir = baseDir.resolve('bin')
quantDir = baseDir.resolve('quantification')

// Check required parameters
if (!params.input) {
    exit 1, "Input file not specified"
}

if (!params.genome) {
    exit 1, "Genome file not specified"
}

if (!params.annotation) {
    exit 1, "Annotation file not specified"
}

log.info "BP Pipeline"
log.info "==========="
log.info "Input                     : ${params.input}"
log.info "Genome                    : ${params.genome}"
log.info "Annotation                : ${params.annotation}"
log.info "Max mismatches            : ${params.mismatches}"
log.info "Max hits                  : ${params.hits}"
log.info "Quality offset            : ${params.quality_offset}"
log.info "Read length               : ${params.max_read_length}"
log.info "Read strandedness         : ${params.read_strand}"
log.info "Log level                 : ${params.loglevel}"
log.info "Number of cpus            : ${params.cpus}"
log.info "Paired                    : ${params.paired_end}"
log.info "Flux count elements       : ${params.count_elements}"
log.info "Read group                : ${params.read_group}"
log.info "Bam stats                 : ${params.bam_stats}"
log.info "Flux memeory              : ${params.flux_mem}"
log.info "Temporary folder          : ${params.tmp_dir}"
log.info "Output folder             : ${file(params.outdir)}"
log.info "Steps to be performed     : ${params.steps.replace(',',' ')}"
log.info "Base folder               : ${baseDir}"
log.info "Bin folder                : ${binDir}"
log.info "Quantification folder     : ${quantDir}"

genome_file = file(params.genome)
annotation_file = file(params.annotation)

input_files = Channel
    .fromPath(params.input) 
    .groupBy {
        path -> path.name.split("\\.", 2)[0][0..-3]
    }
    .flatMap ()
    .map {
       [it.key, it.value.sort()[0], it.value.sort()[1]]
    }

result_path = file(params.outdir)

//result = Channel.create()

//result.subscribe {
//    println it
//}

process index {
    input:
    file genome_file
    
    output:
    file "genome_index.gem" into genome_index
      
    script:
    """
    gemtools index -i ${genome_file} -t ${params.cpus} -o genome_index.gem 
    """
}

(genome_index1, genome_index2) = genome_index.into(2)

process t_index {
    input:
    file genome_index from genome_index1
    file annotation_file

    output:
    set file('tx_index.junctions.gem'), file('tx_index.junctions.keys') into tx_index
    
    script:
    """
    gemtools t-index -i ${genome_index} -a ${annotation_file} -m ${params.max_read_length} -t ${params.cpus} -o tx_index
    """    
}


process mapping {
    input:
    set reads_name, file(read1), file(read2) from input_files
    file genome_index from genome_index2.first()
    set file(tx_index), file(tx_keys) from tx_index.first()

    output:
    set reads_name, view, "mapping.bam" into bam

    script:
    view = 'alignment'
    """
    gemtools rna-pipeline -i ${genome_index} -r ${tx_index} -k ${tx_keys} -f ${read1} -t ${params.cpus} -q ${params.quality_offset} -n mapping
    """
}

(bam1, bam2, bam3, bam4) = bam.into(4)

genomeFai = file(params.genome+'.fai')

process bigwig {
    input:
    set reads_name, in_view, file(bam) from bam1
    file genomeFai

    output:
    set reads_name, view, file('*.bw') into bigwig

    script:
    view = 'bigwig'
    command = ''
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (params.read_strand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        mateBit = (params.read_strand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${params.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${params.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam mapping.bam\n"
    }
    
    strand.each( { 
        command += "genomeCoverageBed " 
        command += (it.key != '' ? "-strand ${it.key} " : '') 
        command += "-split -bg -ibam ${bam} > ${reads_name}${it.value}.bedgraph\n"
        command += "bedGraphToBigWig ${reads_name}${it.value}.bedgraph ${genomeFai} ${reads_name}${it.value}.bw\n" 
    } )

    return command
        
}

process contig {
    input:
    set reads_name, in_view, file(bam) from bam2
    file genomeFai

    output:
    set reads_name, view, file('*_contigs.bed') into contig

    script:
    view = 'contig'
    command = ''
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (params.read_strand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        mateBit = (params.read_strand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${params.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${params.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam mapping.bam\n"
    }

    command += "bamflag -in ${bam} -out tmp.bam -m 3\n"
    command += "mv -f tmp.bam mapping.bam\n"
    
    strand.each( { 
        command += "genomeCoverageBed " 
        command += (it.key != '' ? "-strand ${it.key} " : '') 
        command += "-split -bg -ibam ${bam} > ${reads_name}${it.value}.bedgraph\n"
    } )
    
    if (strand.size() == 2) {
        command += "contigsNew.py --chrFile ${genomeFai}"
        strand.each( {
            command += " --file${it.value.substring(1,2).toUpperCase()} ${reads_name}${it.value}.bedgraph"
        } )
        command += " | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\t\")} print s}'" 
        command += " > ${reads_name}_contigs.bed"
    } else {        
        command += "bamToBed -i ${bam} | sort -k1,1 -k2,2n"
        command += " | mergeBed"
        command += " > ${reads_name}_contigs.bed"
    }
    
    return command

}

process quantification {
    input:
    set reads_name, in_view, file(bam) from bam3
    file annotation_file

    output:
    set reads_name, view, file('flux.gtf') into flux

    script:
    view = 'transcript'
    """
    flux-capacitor -i ${bam} -a ${annotation_file} -o flux.gtf -m AUTO --read-strand ${params.read_strand}
    """
}

process store {

    maxForks 1

    input:
    set reads_name, view, file(store_file) from bam4.mix(bigwig, contig, flux)

    script:
    """
    idxtools add path=`readlink -f ${store_file}` id=${reads_name} view=${view} type=${store_file.name.split("\\.", 2)[1]}
    """
}

//[bam4, bigwig, contig, flux].each {
//    it.subscribe { println it }
//}

//# activate python virtualenv
//run ". $BASEDIR/bin/activate" "$ECHO"
//
//
//
//basename=$(basename $input)
//sample=${basename%.fastq*}
//if [[ $paired == "true" ]];then
//    sample=${basename%[_-\.]1*}
//fi
//
//genomeFai="$genome.fai"
//gemIndex="${genome%.fa}.gem"
//tindex="$annotation.junctions.gem"
//tkeys="$annotation.junctions.keys"
//
//annName=`basename $annotation`
//
//hthreads=$((threads/2))
//if [[ $hthreads == 0 ]];then
//    hthreads=1
//fi
//
//## Binaries
//#
//gem2sam="gem-2-sam"
//samtools="samtools"
//trToGn="TrtoGn_RPKM.sh"
//trToEx="TrtoEx_RPKM.sh"
//bamToContigs="bamToContigs.sh"
//gt_quality="gt.quality"
//gt_filter="gt.filter"
//gt_stats="gt.stats"
//pigz="pigz"
//bamflag="bamflag"
//makecontig="contigsNew.py"
//gtfToGenePred="gtfToGenePred"
//genePredToBed12="genePredToBed12.awk"
//
//## Output files
//#
//# mapping
//gem="$outdir/$sample.map.gz"
//filteredGem=${gem%.map.gz}_m${mism}_n${hits}.map.gz    
//filteredGemStats=${filteredGem%.map.gz}.stats
//filteredBam=${filteredGem%.map.gz}.bam
//filteredBai="$filteredBam.bai"
//# bigwig
//
//# contigs
//contigFile=$outdir/${sample}_contigs.bed
//# quantification
//fluxGtf="$quantDir/$sample.gtf"
//
//## Print pipeline configuration
//
//header="Pipeline configuration for $sample"
//echo $header
//eval "for i in {1..${#header}};do printf \"-\";done"
//printf "\n\n"
//printf "  %-34s %s\n" "Input file:" "$input"
//printf "  %-34s %s\n" "Reference genome file:" "$genome"
//printf "  %-34s %s\n" "Reference gene annotation file:" "$annotation"
//printf "  %-34s %s\n" "Paired-end:" "$paired"
//printf "  %-34s %s\n" "Max number of allowed mismatches:" "$mism"
//printf "  %-34s %s\n" "Max number of hits:" "$hits"
//printf "  %-34s %s\n" "Quality offset:" "$qualityOffset"
//printf "  %-34s %s\n" "Max read length:" "$maxReadLength"
//printf "  %-34s %s\n" "Strandedness:" "$readStrand"
//printf "  %-34s %s\n" "Number of threads:" "$threads"
//printf "  %-34s %s\n" "Flux Capacitor memory:" "$fluxMem"
//printf "  %-34s %s\n" "Temporary folder:" "$tmpdir"
//printf "  %-34s %s\n" "Loglevel:" "$loglevel"
//printf "\n\n"
//printf "  %-34s %s\n" "Mapping folder:" "$outdir"
//printf "  %-34s %s\n" "Quantification folder:" "$quantDir"
//printf "\n\n"
//
//## START
//#
//
//printHeader "Starting Blueprint pipeline for $sample"
//pipelineStart=$(date +%s)
//
//## Mapping
//#
//if [[ $doMapping == "true" ]];then 
//    if [ ! -e $gem ];then
//        step="MAP"
//        startTime=$(date +%s)
//        printHeader "Executing mapping step"                
//            
//        if [ -d $tmpdir ]; then
//            ## Copy needed files to TMPDIR
//            copyToTmp "$gemIndex,$annotation,$tindex,$tkeys"
//            IFS=',' read index annotation tindex tkeys <<< "$paths"
//            gem="$tmpdir/$sample.map.gz"
//        fi
//    
//        log "Running gemtools rna pipeline on ${sample}" $step
//        command="gemtools --loglevel $loglevel rna-pipeline -f $input -q $qualityOffset -i $gemIndex -a $annotation -o `dirname $gem` -t $threads --max-read-length $maxReadLength --no-stats --no-bam"
//        if [[ $paired == "false" ]]; then
//            command="$command --single-end"
//        fi
//        run "$command" "$ECHO"
//   
//        set -e && finalizeStep $gem $tmpdir $outdir
//        IFS=',' read gem <<< "$paths"
//        
//        endTime=$(date +%s)
//        printHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "Map file already present...skipping mapping step"
//    fi
// 
//    ## Filtering the map file
//    ##
//    
//    if [ ! -e $filteredGem ];then
//        step="FILTER"
//        startTime=$(date +%s)
//        printHeader "Executing filtering step"
//        
//        log "Filtering map file..." $step
//        #run "$gt_quality -i $gem -t $threads | $gt_filter --max-levenshtein-error $mism -t $threads | $gt_filter --max-matches $hits -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
//        run "$gt_quality -i $gem -t $threads > $filteredGem.quality.map" "$ECHO"
//        run "$gt_filter -i $filteredGem.quality.map --max-levenshtein-error $mism -t $threads > $filteredGem.levenshtein.map && rm $filteredGem.quality.map" "$ECHO"
//        run "$gt_filter -i $filteredGem.levenshtein.map --max-matches $hits -t $threads | $pigz -p $threads -c > $filteredGem && rm $filteredGem.levenshtein.map" "$ECHO"        
//        log "done\n" $step
//
//        set -e && finalizeStep $filteredGem "-" $outdir
//
//        endTime=$(date +%s)
//        printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "Filtered map file is present...skipping fltering step"
//    fi
// 
//    ## Getting stats for the filtered map file
//    ##
//    
//    if [ $filteredGemStats -ot $filteredGem ]; then
//        step="GEM-STATS"
//        startTime=$(date +%s)
//        printHeader "Executing GEM stats step"
//    
//        log "Producing stats for $filteredGem..." $step
//        command="$gt_stats -i $filteredGem -t $threads -a"
//        if [[ $paired == "true" ]]; then
//            command=$command" -p"
//        fi
//        run "$command 2> $filteredGemStats" "$ECHO"
//        log "done\n" $step
//
//        set -e && finalizeStep $filteredGemStats "-" $outdir
//
//        endTime=$(date +%s)
//        printHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "GEM stats file is present...skipping GEM stats step"
//    fi
//
//    ## Convert to bam 
//    ##
//    
//    if [ ! -e $filteredBam ]; then
//        step="CONVERT"
//        startTime=$(date +%s)
//        printHeader "Executing conversion step"
//
//        if [ -d $tmpdir ]; then
//            ## Copy needed files to TMPDIR
//            copyToTmp "$gemIndex"
//            IFS=',' read index <<< "$paths"
//            filteredBam="$tmpdir/`basename $filteredBam`"
//        fi
//       
//        log "Converting  $sample to bam..." $step
//        command="$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $gemIndex -q offset-$qualityOffset -l"
//        if [[ $readGroup ]]; then
//            command="$command --read-group $readGroup"
//        fi
//        if [[ $paired == "true" ]]; then
//            command="$command --expect-paired-end-reads"
//        else
//            command="$command --expect-single-end-reads | awk 'BEGIN{OFS=FS=\"\t\"}\$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and(\$2,a[i])>0){\$2=xor(\$2,a[i])}}}{print}'"
//        fi
//
//        if [[ $readGroup ]];
//        then
//            command="$command --read-group $readGroup"
//        fi
//        run "$command | $samtools view -@ $threads -Sb - | $samtools sort -@ $threads -m 4G - ${filteredBam%.bam}" "$ECHO"
//        log "done\n" $step
//        
//        set -e && finalizeStep $filteredBam $tmpdir $outdir       
//        IFS=',' read filteredBam <<< "$paths"
//
//        endTime=$(date +%s)
//        printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "Bam file is present...skipping conversion step"
//    fi
//    
//    ## Indexing the filtered bam file
//    ##
//    if [ $filteredBai -ot $filteredBam ];then
//        step="INDEX"
//        startTime=$(date +%s)
//        printHeader "Executing indexing step on the filtered bam file"
//    
//        if [ -d $tmpdir ]; then
//            ## Copy needed files to TMPDIR
//            copyToTmp "$filteredBam"
//            IFS=',' read filteredBam <<< "$paths"
//            filteredBai="$tmpdir/`basename $filteredBai`"
//        fi
//    
//        log "Indexing the filtered bam file\n" $step
//        run "$samtools index $filteredBam" "$ECHO"
//        
//        set -e && finalizeStep $filteredBai $tmpdir $outdir        
//        IFS=',' read filteredBai <<< "$paths"
//    
//        endTime=$(date +%s)
//        printHeader "Indexing step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "Bam index file is present...skipping indexing step"
//    fi
//
//else
//    printHeader "Skipped mapping steps"
//fi
//
//
//## Producing stats for bam file
//##
//statsDir=$outdir/stats
//if [[ $doBamstats == "true" ]]; then
//    step="BAM-STATS"
//    startTime=$(date +%s)
//    printHeader "Executing bam stats step on the filtered bam file"
//
//    if [ ! -d $statsDir ];then
//        log "Creating stats folder..."
//        run "mkdir $statsDir" "$ECHO"
//        log "done\n"
//    fi
//
//    if [ -d $tmpdir ]; then
//        ## Copy needed files to TMPDIR
//        copyToTmp "$filteredBam,$annotation"
//        IFS=',' read filteredBam annotation <<< "$paths"
//    fi
//
//    log "Producing mapping stats for the bam file\n" $step
//
//    ## Create another bam with only uniquely mapping reads
//    uniqBam=${filteredBam%.bam}_uniq.bam
//    if [ ! -e $uniqBam ];then
//        log "Making a bam file of unique mappings..." $step
//        set +e && run "$bamflag -in $filteredBam -out $uniqBam -m 3" "$ECHO"
//        log "done\n"
//    fi
//
//    ## Create a bed12 from gtf
//    genePred=${annotation%.gtf}.GenePred
//    bed12=${annotation%.gtf}.bed12
//    if [ ! -e $genePred ];then
//        set +e && run "$gtfToGenePred $annotation -allErrors $genePred 2> $genePred.err" "$ECHO"
//        set +e && run "cat $genePred | $genePredToBed12 > $bed12" "$ECHO"
//    fi
//
//    ## Check that the files were created successfully
//    if [ ! -e $uniqBam ] || [ ! -e $genePred ];then
//        if [[ ! $ECHO ]];then
//            log "Error producing input file for stats" "ERROR" >&2
//            exit -1
//        fi
//    fi
//
//    ## Run the statistics
//    #
//
//    ## bam_stat
//    bamStat=$statsDir/$sample.bam_stat
//    if [ ! -e $bamStat ];then
//        log "Computing general stats..." $step
//        run "bam_stat.py -i $filteredBam 2> $bamStat" "$ECHO"
//        log "done\n"
//    fi
//
//
//    ## Clipping_profile
//    if [ `ls $statsDir/$sample.clipping_profile.* 2> /dev/null |wc -l` != 3 ];then
//        log "Gettig clipping profile..." $step
//        run "clipping_profile.py -i $filteredBam -o $statsDir/$sample" "$ECHO"
//        log "done\n"
//    fi
//
//    # Gene body coverage
//    # It outputs a file with counts and an R script. If the pdf is not produced,
//    # the R script must be executed manually.
//    # Stderr will contain a line like this:
//    # Cannot generate pdf file from ERR180942_filtered_cuff.geneBodyCoverage_plot.r
//    #
//    if [ `ls $statsDir/$sample.geneBodyCoverage* 2> /dev/null |wc -l` != 3 ];then
//        log "Gettig gene body coverage..." $step
//        run "geneBody_coverage.py -r $bed12 -i $filteredBam -o $statsDir/$sample" "$ECHO"
//        log "done\n"
//    fi
//
//    # Infer the configuration of the experiment, output is on stdout
//    if [ ! -e $statsDir/$sample.inferred_exp ];then
//        log "Inferring the configuration of the experiment..." $step
//        run "infer_experiment.py -r $bed12 -i $filteredBam > $statsDir/$sample.inferred_exp" "$ECHO"
//        log "done\n"
//    fi
//
//    # Splice junctions
//    if [ `ls $statsDir/$sample.[sj]* 2> /dev/null |wc -l` != 5 ];then
//        log "Performing junction annotation..." $step
//        run "junction_annotation.py -i $filteredBam -o $statsDir/$sample -r $bed12 &> $statsDir/$sample.junction_annotation.log" "$ECHO"
//        log "done\n"
//    fi
//
//    if [[ $paired == "true" ]]; then
//        # Insert size distribution
//        if [ `ls $statsDir/$sample.inner_distance* 2> /dev/null |wc -l` != 4 ];then
//            log "Getting the insert size distribution..." $step
//            run "inner_distance.py -i $filteredBam -o $statsDir/$sample -r $bed12" "$ECHO"
//            log "done\n"
//        fi
//    fi
//
//    endTime=$(date +%s)
//    printHeader "BAM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//else
//    printHeader "Skipping BAM stats step"
//fi
//
//
//## Producing bigWig files
//##
//if [[ $doBigWig == "true" ]];then
//    step="BIGWIG"
//    startTime=$(date +%s)
//    printHeader "Executing BigWig step"
//
//    ## Copy needed files to TMPDIR
//    if [ -d $tmpdir ]; then
//        copyToTmp "$filteredBam"
//        IFS=',' read filteredBam <<< "$paths"
//    fi
//
//    log "Producing bigWig files\n" $step
//
//    if [[ $readStrand != "NONE" ]];then
//
//        ## Producing temporary bam with mate1 reversed
//        revBam="${filteredBam%.bam}_1rev.bam"
//        log "Making temporary bam file with mate1 strand reversed..." $step
//        run "$samtools view -h -@ $hthreads $filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $hthreads -Sb - > $revBam" "$ECHO"
//        log "done\n"
//
//        for strand in + -;
//        do
//            suffix="plusRaw"
//            if [[ $strand == "-" ]];then
//                suffix="minusRaw"
//            fi
//
//            bedGraph=$outdir/$sample.$suffix.bedgraph
//            bigWig=$outdir/$sample.$suffix.bw
//            if [ -d $tmpdir ]; then                
//                bedGraph=$tmpdir/$sample.$suffix.bedgraph
//                bigWig=$tmpdir/$sample.$suffix.bw
//            fi
//            log "Making bedGraph $strand strand\n" "$step"
//            run "genomeCoverageBed -strand $strand -split -bg -ibam $revBam > $bedGraph" "$ECHO"
//            log "Making bigWig $strand strand\n" $step
//            run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"
//
//            set -e && finalizeStep $bigWig $tmpdir $outdir
//            IFS=',' read bigWig <<< "$paths"
//        done
//    else
//        bedGraph=$outdir/$sample.bedgraph
//        bigWig=$outdir/$sample.bw
//        if [ -d $tmpdir ]; then                
//            bedGraph=$tmpdir/$sample.bedgraph
//            bigWig=$tmpdir/$sample.bw
//        fi
//
//        log "Making bedGraph\n" "BEDGRAPH"
//        run "genomeCoverageBed -split -bg -ibam $filteredBam > $bedGraph" "$ECHO"
//        log "Making bigWig\n" $step
//        run "bedGraphToBigWig $bedGraph $genomeFai $bigWig" "$ECHO"
//
//        set -e && finalizeStep $bigWig $tmpdir $outdir
//        IFS=',' read bigWig <<< "$paths"
//    fi
//
//    endTime=$(date +%s)
//    printHeader "BigWig step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//else
//    printHeader "Skipping bigwig step"
//fi
//
//## Producing contig files
//##
//if [[ $doContig == "true" ]];then
//    step="CONTIGS"
//    startTime=$(date +%s)
//    printHeader "Executing contigs step"
//
//    if [ -d $tmpdir ]; then
//        ## Copy needed files to TMPDIR
//        copyToTmp "$filteredBam"
//        IFS=',' read filteredBam <<< "$paths"
//        contigFile=$tmpdir/${sample}_contigs.bed
//    fi
//
//    log "Producing contigs file\n" $step
//
//    if [[ $readStrand != "NONE" ]];then
//
//        revBam=${filteredBam%.bam}_1rev.bam
//        if [ ! -e $revBam ];then
//            ## Producing temporary bam with mate1 reversed
//            log "Making temporary bam file with mate1 strand reversed..." $step
//            run "$samtools view -h -@ $threads $filteredBam | awk -v MateBit=64 'BEGIN {OFS=\"\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}' | $samtools view -@ $threads -Sb - > $revBam" "$ECHO"
//            log "done\n"
//        fi
//
//        uniqBam=${revBam%.bam}_uniq.bam
//        if [ ! -e $uniqBam ];then
//            log "Making a bam file of unique mappings..." $step
//            set +e && run "$bamflag -in $revBam -out $uniqBam -m 3" "$ECHO"
//            log "done\n"
//        fi
//
//
//        for strand in + -;
//        do
//            suffix="plusRaw"
//            if [[ $strand == "-" ]];then
//                suffix="minusRaw"
//            fi
//            bedGraph=${uniqBam%.bam}.$suffix.bedgraph
//            log "Making bedGraph $strand strand\n" "$step"
//            run "genomeCoverageBed -strand $strand -split -bg -ibam $revBam > $bedGraph" "$ECHO"
//        done
//
//        log "Generationg the contigs file..." $step
//        run "$makecontig --chrFile $genomeFai --fileP ${uniqBam%.bam}.plusRaw.bedgraph --fileM ${uniqBam%.bam}.minusRaw.bedgraph | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\t\")} print s}' > $contigFile" "$ECHO"
//        log "done\n"
//    else
//        
//        uniqBam=${filteredBam%.bam}_uniq.bam
//        if [ ! -e $uniqBam ];then
//            log "Making a bam file of unique mappings..." $step
//            set +e && run "$bamflag -in $filteredBam -out $uniqBam -m 3" "$ECHO"
//            log "done\n"
//        fi
//        
//        log "Generationg the contigs file..." $step
//        run "bamToBed -i $uniqBam | sort -k1,1 -k2,2n | mergeBed > $contigFile" "$ECHO"
//        log "done\n"
//    fi
//
//    set -e && finalizeStep $contigFile $tmpdir $outdir
//    IFS=',' read contigFile <<< "$paths"
//
//    endTime=$(date +%s)
//    printHeader "Contigs step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//else
//    printHeader "Skipping contigs step"
//fi
//
//## Run FLux
//##
//
//## Run transcript quantification
//#
//if [[ $doFlux == "true" ]];then
//    step="FLUX"
//    
//    export FLUX_MEM=$fluxMem
//    
//    if [ ! -d $quantDir ]; then
//        log "Creating sample folder in $quantDir..." $step
//        run "mkdir -p $quantDir" "$ECHO"
//        log "done\n"
//    fi
//    
//    paramFile="$quantDir/$sample.par"
//    proFile="$quantDir/$sample.profile"
//    
//    # prepare parameter file
//    #
//    # READ_STRAND MATE2_SENSE
//    # ANNOTATION_MAPPING PAIRED_STRANDED
//    # COUNT_ELEMENTS [SPLICE_JUNCTIONS, INTRONS]
//    
//    if [ ! -e $paramFile ]; then
//        run "echo \"# Flux Capacitor parameter file for $sample\" >> $paramFile" "$ECHO"
//        annotationMapping="AUTO"
//        if [[ $readStrand != "NONE" ]]; then
//            run "echo \"READ_STRAND $readStrand\" >> $paramFile" "$ECHO"
//            annotationMapping="STRANDED"
//            if [[ $paired == "true" ]];then 
//                annotationMapping="PAIRED_${annotationMapping}"
//            else
//                annotationMapping="SINGLE_${annotationMapping}"
//            fi
//        fi
//        
//        run "echo \"ANNOTATION_MAPPING $annotationMapping\" >> $paramFile" "$ECHO"    
//        run "echo \"COUNT_ELEMENTS [$countElements]\" >> $paramFile" "$ECHO"
//    fi
//    
//    ## Show Flux parameter file
//    #
//    run "echo \"\"" "$ECHO"
//    run "cat $paramFile" "$ECHO"
//    run "echo \"\"" "$ECHO"
//
//    startTime=$(date +%s)
//    printHeader "Executing Flux quantification step"
//
//    if [ -d $tmpdir ];then
//        ## Copy needed files to TMPDIR
//        copyToTmp "$filteredBam,$filteredBai,$annotation"
//        IFS=',' read filteredBam filteredBai annotation <<< "$paths"
//        fluxGtf="$tmpdir/$sample.gtf"
//    fi
//
//    if [ ! -e ${annotation%.gtf}_sorted.gtf ];then
//        sortLog="$quantDir/${sample}_sort_annotation.log"
//        log "Checking if the annotation is sorted" $step
//        set -e && run "flux-capacitor -t sortGTF -c -i $annotation -o ${annotation%.gtf}_sorted.gtf > $sortLog 2>&1" "$ECHO"
//        log "done\n"
//    fi
//
//    anno=$annotation
//    if [ -e {$anno%.gtf}_sorted.gtf ];then
//        anno=${anno%.gtf}_sorted.gtf
//        log "Using sorted annotation: $anno\n" $step
//    fi
//
//    if [[ ! -e $proFile ]];then
//        profileLog="$quantDir/${sample}_flux_profile.log"
//        log "Getting sistematic biases along transcripts\n" $step
//        run "flux-capacitor --profile -p $paramFile -i $filteredBam -a $anno --profile-file $proFile > $profileLog 2>&1" "$ECHO"
//    fi
//
//    if [ -e $fluxGtf ];then
//        log "Removing old quantification file" $step
//        rm $fluxGtf
//    fi
//
//    log "Running Flux Capacitor\n" $step
//    quantLog=$quantDir/${sample}_flux_quantification.log
//    run "flux-capacitor -p $paramFile -i $filteredBam -a $anno -o $fluxGtf > $quantLog 2>&1" "$ECHO"
//
//    set -e && finalizeStep $fluxGtf $tmpdir $quantDir
//    IFS=',' read fluxGtf <<< "$paths"
//
//    endTime=$(date +%s)
//    printHeader "Quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//
//    txFile=$quantDir/${sample}_transcript.gtf
//    if [ ! -e $txFile ];then
//        startTime=$(date +%s)
//        printHeader "Getting transcript from quantifications"
//        log "Generating transcripts file..."
//        run "awk '\$3==\"transcript\"' $fluxGtf > $txFile" "$ECHO"
//        log "done\n"
//        log "Computing md5sum for transcripts file..." $step
//        run "md5sum $txFile > $txFile.md5" "$ECHO"
//        log "done\n"
//        endTime=$(date +%s)
//        printHeader "Transcripts written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    else
//        printHeader "Transcript file present...skipping transcript step"
//    fi
//    
//    if [[ $countElements =~ "SPLICE_JUNCTIONS" ]]; then
//        junctionFile=$quantDir/${sample}_junction.gtf
//        if [ ! -e $junctionFile ];then
//            startTime=$(date +%s)
//            printHeader "Getting junctions from quantifications"
//            log "Generating junctions file..."
//            run "awk '\$3==\"junction\"' $fluxGtf > $junctionFile" "$ECHO"
//            log "done\n"
//            log "Computing md5sum for junctions file..." $step
//            run "md5sum $junctionFile > $junctionFile.md5" "$ECHO"
//            log "done\n"
//            endTime=$(date +%s)
//            printHeader "Junctions written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//        else
//            printHeader "Junctions file present...skipping transcript step"
//        fi
//    fi
//    
//    if [[ $countElements =~ "INTRONS" ]]; then
//        intronFile=$quantDir/${sample}_intron.gtf
//        if [ ! -e $intronFile ];then
//            printHeader "Getting all-intronic regions from quantifications"
//            log "Generating introns file..."
//            run "awk '\$3==\"intron\"' $fluxGtf > $intronFile" "$ECHO"
//            log "done\n"
//            log "Computing md5sum for introns file..." $step
//            run "md5sum $intronFile > $intronFile.md5" "$ECHO"
//            log "done\n"
//            endTime=$(date +%s)
//            printHeader "Introns written in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//        else
//            printHeader "Introns file present...skipping transcript step"
//        fi
//    fi
//    
//    exonFile=$quantDir/${sample}_exon_distinct_with_rpkm.gff
//    if [ ! -e $exonFile ];then
//        step="EXON"
//        startTime=$(date +%s)
//        printHeader "Executing Exon quantification step"
//    
//        
//        if [ -d $tmpdir ]; then
//            ## Copy needed files to TMPDIR
//            copyToTmp "$annotation,$fluxGtf"
//            IFS=',' read annotation fluxGtf <<< "$paths"
//            exonFile=$tmpdir/${sample}_exon_distinct_with_rpkm.gff
//        fi
//    
//        log "Running Exon quantification\n" $step
//        run "$trToEx -a $annotation -i $fluxGtf -o `dirname $exonFile`" "$ECHO"
//    
//        set -e && finalizeStep $exonFile $tmpdir $quantDir
//        IFS=',' read exonFile <<< "$paths"
//        
//        endTime=$(date +%s)
//        printHeader "Exon quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    
//    else
//        printHeader "Exon quantification file present...skipping Exon quantification step"
//    fi
//    
//    geneFile=$quantDir/${sample}_gene_with_rpkm.gff
//    if [ ! -e $geneFile ];then
//        step="GENE"
//        startTime=$(date +%s)
//        printHeader "Executing Gene quantification step"
//    
//        if [ -d $tmpdir ]; then
//            ## Copy needed files to TMPDIR
//            copyToTmp "$annotation,$fluxGtf"
//            IFS=',' read annotation fluxGtf <<< "$paths"
//            geneFile=$tmpdir/${sample}_gene_with_rpkm.gff
//        fi
//    
//        log "Running Gene quantification\n" $step
//        run "$trToGn -a $annotation -i $fluxGtf -o `dirname $geneFile`" "$ECHO"
//    
//        set -e && finalizeStep $geneFile $tmpdir $quantDir
//        IFS=',' read geneFile <<< "$paths"
//        
//        endTime=$(date +%s)
//        printHeader "Gene quantificaton step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
//    
//    else
//        printHeader "Gene quantification file present...skipping Gene quantification step"
//    fi
//else
//    printHeader "Skipping Flux quantification step"
//fi
//
//
//# deactivate python virtualenv
//run "deactivate" "$ECHO"
//
//pipelineEnd=$(date +%s)
//
//log "\n"
//printHeader "Blueprint pipeline for $sample completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
//
//# disable extglob
//shopt -u extglob
//
//exit 0
