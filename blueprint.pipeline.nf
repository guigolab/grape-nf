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
params.maxMismatches = 4
params.maxMultimaps = 10
params.bamStats = false
params.countElements = [] 
params.fluxMem = '3G'
params.fluxProfile = false 
if ( params.chunkSize != null ) params.chunkSize = params.chunkSize as int

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

//print usage
if (params.help) {
    log.info '''
B L U E P R I N T ~ RNA Pipeline
--------------------------------
Run the RNAseq pipeline on one sample.

Usage: 
    blueprint.pipeline.nf -i INDEX_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...

Options:
    --help                              Show this message and exit.
    --index INDEX_FILE                  Index file.
    --genome GENOME_FILE                Reference genome file(s).
    --annotation ANNOTAION_FILE         Reference gene annotation file(s).
    --tmp-dir                           Specify the temporary folder to be used as a scratch area.
                                        Default: "$TMPDIR" if the environment variable is defined, "-" otherwise.
    --paired-end                        Specify whether the data is paired-end. Default: "auto".
    --max-read-length READ_LENGTH       The maximum read length (used to compute the transcriptomes). Default: "150".
    --max-mismatches THRESHOLD          Set maps with more than <threshold> error events to unmapped. Default "4".
    --max-multimaps THRESHOLD           Set multi-maps with more than <threshold> mappings to unmapped. Default "10".
    --filter-intron-length THRESHOLD    Filter multimaps preferring ones with intron length > <threshold>
    --filter-block-length THRESHOLD     Filter multimaps preferring ones with block length > <threshold>
    --filter-level THRESHOLD            Reduce multimaps using the specified uniqueness level.
    --read-strand READ_STRAND           Directionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default "auto".
    --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.
    --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.
    --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.
    --rg-desc DESCRIPTION               Description for the BAM @RG tag.
    --flux-mem MEMORY                   Specify the amount of ram the Flux Capacitor can use. Default: "3G".
    --flux-profile                      Specify whether the Flux Capacitor profile file shoudl be written. Default: "false".
    --count-elements ELEMENTS           A comma separated list of elements to be counted by the Flux Capacitor.
                                        Possible values: INTRONS,SPLICE_JUNCTIONS. Defalut: "none".
    --loglevel LOGLEVEL                 Log level (error, warn, info, debug). Default "info".
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
log.info "Input chunk size          : ${params.chunkSize != null ? params.chunkSize : 'no split'}"
log.info "Use temporary folder      : ${params.tmpDir}"
log.info ""
if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Max mismatches            : ${params.maxMismatches}"
    log.info "Max multimaps             : ${params.maxMultimaps}"
    log.info "Max read length           : ${params.maxReadLength != null ? params.maxReadLength : 'auto'}"
    log.info "Read strandedness         : ${params.readStrand != null ? params.readStrand : 'auto'}"
    log.info "Paired                    : ${params.pairedEnd != null ? params.pairedEnd : 'auto'}"
    log.info "Produce BAM stats         : ${params.bamStats}"
    if ( params.rgPlatform ) log.info "Sequencing platform       : ${params.rgPlatform}"  
    if ( params.rgLibrary ) log.info "Sequencing library        : ${params.rgLibrary}"  
    if ( params.rgCenterName ) log.info "Sequencing center         : ${params.rgCenterName}" 
    if ( params.rgDesc ) log.info "@RG Descritpiton          : ${params.rgDesc}" 
    log.info ""
}
if ('flux' in pipelineSteps || 'quantification' in pipelineSteps) {
    log.info "Flux Capacitor parameters"
    log.info "-------------------------"
    log.info "Additional quantified elements : ${params.countElements.size()==0 ? 'NONE' : params.countElements.join(" ")}"
    log.info "Memory                         : ${params.fluxMem}"
    log.info "Create profile file            : ${params.fluxProfile}"
    log.info ""
}

genomes=params.genome.split(',').collect { file(it) }
annos=params.annotation.split(',').collect { file(it) }

index = params.index ? file(params.index) : System.in
input_files = Channel.create()
input_chunks = Channel.create()

Channel
    .from(index.readLines())
    .map {
        line -> [ line.split()[0], line.split()[1], file(line.split()[2]), line.split()[3], line.split()[4] ]
    }
    .subscribe onNext: { sample, id, path, type, view -> if( params.chunkSize != null ) { file(path).splitFastq(by: (params.chunkSize), file: true, autoClose: false, into: input_chunks) { chunk, source -> tuple(sample, id+chunk.baseName.find(/\..+$/), chunk, type, view) } } else {input_chunks << tuple(sample, id, path, type, view)} }, onComplete: { input_chunks << Channel.STOP }

input_files = input_chunks
    .groupBy {
        sample, id, path, type, view -> id 
    }
    .flatMap ()
    .map {        
       [it.key, it.value[0][0], it.value.collect { sample, id, path, type, view -> path }, fastq(it.value[0][2]).qualityScore()]
   }

Genomes = Channel.create()
Annotations = Channel.create()
Channel.from(genomes)
    .merge(Channel.from(annos)) {
        g, a -> [g.name.split("\\.", 2)[0], g, a]
    }
    .separate(Genomes, Annotations) {
        a -> [[a[0],a[1]], [a[0],a[2]]]
    }

(Genomes1, Genomes2) = Genomes.into(2)
(Annotations1, Annotations2, Annotations3, Annotations4, Annotations5, Annotations6) = Annotations.into(6)

process fastaIndex {
    input:
    set species, file(genome) from Genomes1
    set species, file(annotation) from Annotations1

    output:
    set species, file("${genome}.fai") into FaiIdx

    script:
    def command = ""
    
    command += "samtools faidx ${genome}"
    
    return command
}

(FaiIdx1, FaiIdx2) = FaiIdx.into(2)

process index {

    input:
    set species, file(genome) from Genomes2

    output:
    set species, file("genome_index.gem") into GenomeIdx

    script:
    def command = ""
    
    command += "gemtools index -i ${genome} -t ${task.cpus} -o genome_index.gem"
    
    return command
}

(GenomeIdx1, GenomeIdx2, GenomeIdx3) = GenomeIdx.into(3)

process t_index {

    input:
    set species, file(genome_index) from GenomeIdx1
    set species, file(annotation) from Annotations2

    output:
    set species, file('tx_index.junctions.gem'), file('tx_index.junctions.keys') into TranscriptIdx

    script:
    def command = ""

    command += "gemtools t-index -i ${genome_index} -a ${annotation} -t ${task.cpus} -o tx_index" 
    if ( params.maxReadLength != null ) command += ' -m ${params.maxReadLength}'

    return command
}

pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

process mapping {

    input:
    set id, sample, file(reads), qualityOffset from input_files
    set species, file(annotation) from Annotations3.first()
    set species, file(genome_index) from GenomeIdx2.first()
    set species, file(tx_index), file(tx_keys) from TranscriptIdx.first()

    output:
    set id, sample, view, "${id}.filtered.map.gz", pairedEnd, qualityOffset into map

    script:
    view = 'gemFiltered'
    def command = ""
    
    fqs = reads.toString().split(" ")
    pairedEnd = false
    if (fqs.size() == 2) pairedEnd = true 

    command += "gemtools rna-pipeline -i ${genome_index} -a ${annotation} -r ${tx_index} -k ${tx_keys} -f ${reads}"
    command += " --filter-max-multi-maps ${params.maxMultimaps}"
    command += " --filter-max-error-events ${params.maxMismatches}"
    if ( params.filterBlockLength ) command += " --filter-block-length ${params.filterBlockLength}"
    if ( params.filterIntronLength ) command += " --filter-intron-length ${params.filterIntronLength}"
    if ( params.filterUniqLevel ) command += " --filter-level ${params.filterUniqLevel}"
    command += " --no-bam"
    if (!pairedEnd) {
        command += " --single-end"
    }
    command += " -t ${task.cpus} -q ${qualityOffset} -n ${id}"

    return command
}

process gemToBam {

    input:
    set id, sample, view, gem_filtered, pairedEnd, qualityOffset from map
    set species, file(genome_index) from GenomeIdx3.first()

    output:
    set id, sample, type, view, "${id}${prefix}.bam", pairedEnd into bam

    script:    
    // prepare BAM @RG tag information
    // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
    def date = ""
    def readGroup = []
    readGroup << "ID=${id}" 
    readGroup << "PU=${id}" 
    readGroup << "SM=${sample}" 
    if ( date ) readGroup << "DT=${date}"
    if ( params.rgPlatform ) readGroup << "PL=${params.rgPlatform}"
    if ( params.rgLibrary ) readGroup << "LB=${params.rgLibrary}"
    if ( params.rgCenterName ) readGroup << "CN=${params.rgCenterName}"
    if ( params.rgDesc ) readGroup << "DS=${params.rgDesc}"

    def command = ""
    awkCommand = 'BEGIN{OFS=FS=\"\t\"}\$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and(\$2,a[i])>0){\$2=xor(\$2,a[i])}}}{print}'
    type = "bam"
    view = "Alignments"
    prefix = pref

    command += "pigz -p ${task.cpus} -dc ${gem_filtered}"
    command += " | gem-2-sam -T ${task.cpus} -I ${genome_index} -q offset-${qualityOffset} -l"
    if (readGroup) {
       command += " --read-group ${readGroup.join(',')}"
    }
    if (pairedEnd) {
       command += " --expect-paired-end-reads"
    }
    else {
       command += " --expect-single-end-reads"
       command += " | awk '${awkCommand}'"
    }

    command += " | samtools view -@ ${task.cpus} -Sb -"
    command += " | samtools sort -@ ${task.cpus} -m ${(long)(task.memory.toBytes()/(2*task.cpus))} - ${id}${prefix}"
    command += " && samtools index ${id}${prefix}.bam"

    return command
}

singleBam = Channel.create()
groupedBam = Channel.create()

bam.groupBy() {
    id, sample, type, view, path, pairedEnd -> sample
}
.flatMap()
.map {
    [it.value.size() > 1 ? it.key : it.value[0][0], it.value[0][2], it.value[0][3], it.value.collect { id, sample, type, view, path, pairedEnd -> path }, it.value[0][5]]
}.choice(singleBam, groupedBam) {
    it[3].size() > 1 ? 1 : 0
}

process mergeBam {
    
    input:
    set id, type, view, file(bam), pairedEnd from groupedBam

    output:
    set id, type, view, "${id}${prefix}.bam", pairedEnd into mergedBam

    script:
    def command = ""
    prefix = pref

    command += "(samtools view -H ${bam} | grep -v \"@RG\";for f in ${bam};do samtools view -H \$f | grep \"@RG\";done) > header.txt\n"
    command += "samtools merge -@ ${task.cpus} -h header.txt ${id}${prefix}.bam ${bam}"
}

bam = singleBam
.mix(mergedBam)
.map { 
    it.flatten() 
}

(bam1, bam2) = bam.into(2)

process inferExp {
    input:
    set id, type, view, file(bam), pairedEnd from bam1
    set species, file(annotation) from Annotations4.first()

    output:
    set id, stdout into bamInf

    script:
    def command = ""
    def genePred = "${annotation.name.split('\\.', 2)[0]}.genePred"
    def bed12 = "${annotation.name.split('\\.', 2)[0]}.bed"

    command += "set -o pipefail\n"
    command += "gtfToGenePred ${annotation} -allErrors -ignoreGroupsWithoutExons ${genePred} 2> ${genePred}.err\n"
    command += "genePredToBed ${genePred} ${bed12}\n" 
    command += "${baseDir}/bin/infer_experiment.py -i ${bam} -r ${bed12} 2> infer_experiment.log | tr -d '\\n'"
}

(bamInf1, bamInf2) = bamInf.into(2)

bamInf1.subscribe {        
    tuple ->
        (id, readStrand) = tuple
        if (params.readStrand != null && !readStrand.equals(params.readStrand)) {
            log.warn "----> '${id}' skipped"
            log.warn "Detected and supplied read strandedness do not match:"
            log.warn "${readStrand} != ${params.readStrand}"
        }
}

bamStrand = bam2.phase(bamInf2)
.map {
    [it[0],it[1][1]].flatten()
}.filter { id, type, view, bam, pairedEnd, readStrand ->
    params.readStrand == null || readStrand.equals(params.readStrand)
}

(bam1, bam2, bam3, out) = bamStrand.into(4)

process bigwig {
    
    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam1
    set species, file(genomefai) from FaiIdx1.first()
    
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
        command += "bedGraphToBigWig ${id}${it.value}.bedgraph ${genomefai} ${id}${it.value}.bw\n"

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
    set species, file(genomefai) from FaiIdx2.first()

    output:
    set id, type, view, file('*_contigs.bed') into contig

    script:
    type = 'bed'
    view = 'Contigs'
    bamflagMode = '2'
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

    command += "bamflag -in ${bam} -out tmp.bam -m ${bamflagMode}\n"
    command += "mv -f tmp.bam ${bam}\n"

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${id}${it.value}.bedgraph\n"
    } )

    if (strand.size() == 2) {
        command += "contigsNew.py --chrFile ${genomefai}"
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
    set id, type, view, file(bam), file(bai), pairedEnd, readStrand from bam3.map { [it[0],it[1],it[2],it[3],file("${it[3].toAbsolutePath()}.bai"),it[4],it[5]] }
    set species, file(annotation) from Annotations5.first() 

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

(flux1, flux2) = flux.into(2)

process geneQuantification {
    
    input:
    set id, type, view, file(fluxGtf) from flux1
    set species, file(annotation) from Annotations6.first() 

    output:
    set id, type, view, file("${id}${prefix}_gene_with_rpkm.gff") into genes

    script:
    type = "gff"
    view = "Gene${annotation.name.replace('.gtf','').capitalize()}"
    prefix = pref
    command = "TrtoGn_RPKM.sh -a ${annotation} -i ${fluxGtf}"
    
    return command

}

out.mix(bigwig, contig, flux2, genes).collectFile(name: "pipeline.db", newLine: true) {
    [it[3], it[0], it[1], it[2]].join("\t")
}
.subscribe {
    def msg = "Output files db -> ${it}"
    log.info ""
    log.info "-" * msg.size()
    log.info "Pipeline run completed."
    log.info ""
    log.info msg
    log.info "-" * msg.size()
}
