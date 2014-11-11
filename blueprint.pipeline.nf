#!/bin/env nextflow 
/*
 * Copyright (c) 2014, Centre for Genomic Regulation (CRG)
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

//Set default values for params

params.steps = 'mapping,bigwig,contig,flux'
params.tmpDir = (System.env.TMPDIR != null ? true : false)
params.mismatches = 4
params.hits = 10
params.qualityOffset = 33
params.maxReadLength = 150
params.bamStats = false
params.countElements = []
params.fluxMem = '3G'
params.fluxProfile = false 

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

//print usage
if (params.help) {
    log.info '''
B L U E P R I N T ~ RNA Pipeline
--------------------------------
Run the RNAseq pipeline on one sample.

Usage: 
    ./blueprint.pipeline.sh -i FASTQ_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...

Options:
    --help              Show this message and exit.
    --index             Index file.
    --genome            Reference genome file(s).
    --annotation        Reference gene annotation file(s).
    --tmp-dir           Specify the temporary folder to be used as a scratch area.
                        Default: "$TMPDIR" if the environment variable is defined, "-" otherwise.
    --paired-end        Specify whether the data is paired-end. Default: "auto".
    --mismatches        Max number of mismatches. Default "4".
    --hits              Max number of hits. Default "10".
    --quality-offset    The quality offset of the fastq files. Default: "33".
    --max-read-length   The maximum read length (used to compute the transcriptomes). Default: "150".
    --read-strand       directionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default "NONE".
    --flux-mem          Specify the amount of ram the Flux Capacitor can use. Default: "3G".
    --flux-profile      Specify whether the Flux Capacitor profile file shoudl be written. Default: "false".
    --count-elements    A comma separated list of elements to be counted by the Flux Capacitor.
                        Possible values: INTRONS,SPLICE_JUNCTIONS. Defalut: "none".
    --loglevel          Log level (error, warn, info, debug). Default "info".
    '''
    exit 1
}

if (!params.genome) {
    exit 1, "Genome file not specified"
}

if (!params.annotation) {
    exit 1, "Annotation file not specified"
}

log.info ""
log.info "B L U E P R I N T ~ RNA Pipeline"
log.info ""
log.info "General parameters"
log.info "------------------"
log.info "Index file                : ${params.index}"
log.info "Genome                    : ${params.genome}"
log.info "Annotation                : ${params.annotation}"
log.info "Steps to be performed     : ${pipelineSteps.join(" ")}"
log.info "Use temporary folder      : ${params.tmpDir}"
log.info ""
if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Max mismatches            : ${params.mismatches}"
    log.info "Max multimaps             : ${params.hits}"
    log.info "Quality offset            : ${params.qualityOffset}"
    log.info "Read length               : ${params.maxReadLength}"
    log.info "Read strandedness         : ${params.readStrand != null ? params.readStrand : 'auto'}"
    log.info "Paired                    : ${params.pairedEnd != null ? params.pairedEnd : 'auto'}"
    log.info "Read group                : ${params.readGroup != null ? params.readGroup : '-'}"
    log.info "Produce BAM stats         : ${params.bamStats}"
    log.info ""
}
if ('flux' in pipelineSteps || 'quantification' in pipelineSteps) {
    log.info "Flux Capacitor parameters"
    log.info "-------------------------"
    log.info "Elements to be quantified : ${params.countElements.join(" ")}"
    log.info "Memory                    : ${params.fluxMem}"
    log.info "Create profile file       : ${params.fluxProfile}"
    log.info ""
}

genomes=params.genome.split(',').collect { file(it) }
genomeFais=params.genome.split(',').collect { file("${it}.fai") }
annos=params.annotation.split(',').collect { file(it) }

index = params.index ? file(params.index) : System.in

input_files = Channel
    .from(index.readLines())
    .map {
        line -> [ file(line.split()[0]), line.split()[1], line.split()[2], line.split()[3] ]
    }
    .groupBy {
        path, id, type, view -> id 
    }
    .flatMap ()
    .map {
       [it.key, it.value.collect { path, id, type, view -> path }]//.flatten()
    }

Refs = Channel.from(genomes)
    .merge(Channel.from(annos)) {
        g, a -> [g.name.split("\\.", 2)[0], g, a]
    }

process index {

    input:
    set species, file(genome), file(annotation) from Refs

    output:
    set species, file(genome), file(annotation), file("genome_index.gem") into Refs

    script:
    view = "genomeIdx"
    def command = ""
    
    command += "gemtools index -i ${genome} -t ${task.cpus} -o genome_index.gem"
    
    return command
}

process t_index {

    input:
    set species, file(genome), file(annotation), file(genome_index) from Refs

    output:
    set species, file(genome), file(annotation), file(genome_index), file('tx_index.junctions.gem'), file('tx_index.junctions.keys') into IdxRefs

    script:
    def command = ""

    command += "gemtools t-index -i ${genome_index} -a ${annotation} -m ${params.maxReadLength} -t ${task.cpus} -o tx_index" 

    return command
}

(IdxRefs1, IdxRefs2) = IdxRefs.into(2)

process mapping {

    input:
    set id, file(reads) from input_files
    set species, file(genome), file(annotation), file(genome_index), file(tx_index), file(tx_keys) from IdxRefs1.first()

    output:
    set id, view, "${id}.map.gz", pairedEnd into map

    script:
    view = 'gemUnfiltered'
    def command = ""
    
    fqs = reads.toString().split(" ")
    pairedEnd = false
    if (fqs.size() == 2) pairedEnd = true 

    command += "gemtools rna-pipeline -i ${genome_index} -a ${annotation} -r ${tx_index} -k ${tx_keys} -f ${fqs[0]}"
    command += " --no-stats --no-bam"
    if (!pairedEnd) {
        command += " --single-end"
    }
    command += " -t ${task.cpus} -q ${params.qualityOffset} -n ${id}"

    return command
}

pref = "_m${params.mismatches}_n${params.hits}"

process filter {

    input:
    set id, view, gem_unfiltered, pairedEnd from map

    output:
    set id, view, file("${id}${prefix}.map.gz"), pairedEnd into fmap

    script:
    view = "gemFiltered"
    prefix = pref
    def command = ""

    command += "${baseDir}/bin/gt.quality -i ${gem_unfiltered} -t ${task.cpus} > quality.map;"
    command += "${baseDir}/bin/gt.filter -i quality.map --max-levenshtein-error ${params.mismatches} -t ${task.cpus} > mism.map && rm quality.map;"
    command += "${baseDir}/bin/gt.filter -i mism.map --max-matches ${params.hits} -t ${task.cpus} > multimaps.map && rm mism.map;"
    command += "pigz -p ${task.cpus} -c multimaps.map > ${id}_m${params.mismatches}_n${params.hits}.map.gz && rm multimaps.map"

    return command
}

(fmap1, fmap2) = fmap.into(2)

process gemStats {

    input:
    set id, view, gem_filtered, pairedEnd from fmap1

    output:
    set id, view, "${id}${prefix}.stats", pairedEnd into stats

    script:
    def command = ""
    view = "gemFilteredStats"
    prefix = pref
    
    command += "gt.stats -i ${gem_filtered} -t ${task.cpus} -a"
    if (pairedEnd) {
        command += " -p"
    }
    command += " 2> ${id}${prefix}.stats"

    return command
}

process gemToBam {

    input:
    set id, view, gem_filtered, pairedEnd from fmap2
    set species, file(genome), file(annotation), file(genome_index), file(tx_index), file(tx_keys) from IdxRefs2.first()

    output:
    set id, type, view, "${id}${prefix}.bam", pairedEnd into bam

    script:
    def command = ""
    awkCommand = 'BEGIN{OFS=FS=\"\t\"}\$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and(\$2,a[i])>0){\$2=xor(\$2,a[i])}}}{print}'
    type = "bam"
    view = "Alignments"
    prefix = pref

    command += "pigz -p ${task.cpus} -dc ${gem_filtered}"
    command += " | gem-2-sam -T ${task.cpus} -I ${genome_index} -q offset-${params.qualityOffset} -l"
    if (params.readGroup) {
       command += " --read-group ${params.readGroup}"
    }
    if (pairedEnd) {
       command += " --expect-paired-end-reads"
    }
    else {
       command += " --expect-single-end-reads"
       command += " | awk '${awkCommand}'"
    }

    command += " | samtools view -@ ${task.cpus} -Sb -"
    command += " | samtools sort -@ ${task.cpus} -m 4G - ${id}${prefix}"
    command += " && samtools index ${id}${prefix}.bam"

    return command
}

process inferExp {
    input:
    set id, type, view, file(bam), pairedEnd from bam
    file anno from Channel.from(annos).first()

    output:
    set id, type, view, file(bam), pairedEnd, stdout into bamInf

    script:
    def command = ""
    def genePred = "${anno.name.split('\\.', 2)[0]}.genePred"
    def bed12 = "${anno.name.split('\\.', 2)[0]}.bed"

    command += "gtfToGenePred ${anno} -allErrors -ignoreGroupsWithoutExons ${genePred} 2> ${genePred}.err;"
    command += "genePredToBed12.awk ${genePred} > ${bed12};" 
    command += "set +eu; . /software/rg/el6.3/virtualenvs/python2.7.3/bin/activate;"
    command += "set -eu; ${baseDir}/bin/infer_experiment.py -i ${bam} -r ${bed12} 2> /dev/null | tr -d '\n'"
}

(bamInf1, bamInf2) = bamInf.into(2)
bamInf1.subscribe {        
    tuple ->
        (id, type, view, bam, pairedEnd, readStrand) = tuple
        if (params.readStrand != null && !readStrand.equals(params.readStrand)) {
            log.warn "----> '${id}' skipped"
            log.warn "Detected and supplied read strandedness do not match:"
            log.warn "${readStrand} != ${params.readStrand}"
        }
}

bamStrand = bamInf2.filter { id, type, view, bam, pairedEnd, readStrand ->
        params.readStrand == null || readStrand.equals(params.readStrand)
}

fais = Channel.from(genomeFais)
(bam1, bam2, bam3, bam4) = bamStrand.into(4)
(fais1, fais2) = fais.into(2)

process bigwig {
    
    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam1
    file genomeFai from fais1.first()

    output:
    set id, type, views, file('*.bw') into bigwig

    script:
    views = []
    view = 'Signal'
    type = "bigwig"
    def command = ''
    strand = ['': '.raw']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (readStrand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        if (pairedEnd) mateBit = (readStrand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${task.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${task.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam ${bam}\n"
    }

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${id}${it.value}.bedgraph\n"
        command += "bedGraphToBigWig ${id}${it.value}.bedgraph ${genomeFai} ${id}${it.value}.bw\n"

        views << "${it.value[1..-1].capitalize()}${view}"
    } )

    return command

}



bigwig = bigwig.reduce([:]) { files, tuple ->
    (id, type, view, path) = tuple     
    if (!files) files = []
    paths = path.toString().replaceAll(/[\[\],]/,"").split(" ")
    (1..paths.size()).each { files << [id, type, view[it-1], paths[it-1]] }
    return files
}
.flatMap()

process contig {

    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam2
    file genomeFai from fais2.first()

    output:
    set id, type, view, file('*_contigs.bed') into contig

    script:
    type = 'bed'
    view = 'Contigs'
    def command = "" 
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (readStrand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        if (pairedEnd) mateBit = (readStrand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${task.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${task.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam ${bam}\n"
    }

    command += "bamflag -in ${bam} -out tmp.bam -m 3\n"
    command += "mv -f tmp.bam ${bam}\n"

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${id}${it.value}.bedgraph\n"
    } )

    if (strand.size() == 2) {
        command += "contigsNew.py --chrFile ${genomeFai}"
        strand.each( {
            command += " --file${it.value.substring(1,2).toUpperCase()} ${id}${it.value}.bedgraph"
        } )
        command += " | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\\t\")} print s}'"
        command += " > ${id}_contigs.bed"
    } else {
        command += "bamToBed -i ${bam} | sort -k1,1 -k2,2n"
        command += " | mergeBed"
        command += " > ${id}_contigs.bed"
    }
    
    return command

}

process quantification {

    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam3
    file annotation from Channel.from(annos).first() 

    output:
    set id, type, view, file("${id}${prefix}.gtf") into flux

    script:
    type = "gtf"
    view = "Transcript${annotation.name.replace('.gtf','').capitalize()}"
    prefix = pref
    def command = ""
    paramFile = file('${id}${prefix}.par')

    //Flux parameter file has to be in the same folder as the bam - Flux bug

    /* paramFile.append("# Flux Capacitor parameter file for ${id}\n")
    annotationMapping = "AUTO"
    if (readStrand != "NONE") {        
        paramFile.append("READ_STRAND ${readStrand}\n")
        annotationMapping="STRANDED"
        if (pairedEnd) {
            annotationMapping="PAIRED_${annotationMapping}"
        }
        else {
            annotationMapping="SINGLE_${annotationMapping}"
        }
    }
    paramFile.append("ANNOTATION_MAPPING ${annotationMapping}\n")
    paramFile.append("COUNT_ELEMENTS ${params.countElements}\n")

    if (params.fluxProfile) {
        paramFile.append("PROFILE_FILE ${id}${prefix}_profile.json\n")
        command += "flux-capacitor --profile -p ${paramFile} -i ${bam}  -a ${annotation}"
    }
    command += "flux-capacitor -p ${paramFile} -i ${bam} -a ${annotation} -o ${id}${prefix}.gtf".toString() */

    // Workaround
    fluxParams = ""
    annotationMapping = "AUTO"
    if (readStrand != "NONE") {
        annotationMapping="STRANDED"
        if (pairedEnd) {
            annotationMapping="PAIRED_${annotationMapping}"            
        }
        else {
            annotationMapping="SINGLE_${annotationMapping}"
        }
        fluxParams += " --read-strand ${readStrand}"
    }
    fluxParams += " -m ${annotationMapping}"
    fluxParams += " --count-elements ${params.countElements}"

    if (params.fluxProfile) {
        fluxParams += " --profile-file ${id}${prefix}_profile.json"
        command += "flux-capacitor --profile ${fluxParams} -i ${bam}  -a ${annotation}; "
    }
    command += "flux-capacitor ${fluxParams} -i ${bam} -a ${annotation} -o ${id}${prefix}.gtf"

    return command
}

bam4.mix(bigwig, contig, flux).collectFile(name: "pipeline.db", newLine: true) {
    [it[3], it[0], it[1], it[2]].join("\t")
}
.subscribe {
    def msg = "Output files db -> ${it}"
    log.info ""
    log.info "-" * msg.size()
    log.info "Pipeline run completed."
    log.info ""
    log.info msg
    //log.info "${it}"
    log.info "-" * msg.size()
}
