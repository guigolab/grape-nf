#!/bin/env nextflow 
/*
 * Copyright (c) 2015, Centre for Genomic Regulation (CRG)
 * Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
 *
 * This file is part of the GRAPE RNAseq pipeline.
 *
 * The GRAPE RNAseq pipeline is a free software: you can redistribute it 
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
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
params.steps = 'mapping,bigwig,contig,quantification'
params.maxMismatches = 4
params.maxMultimaps = 10

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

//print usage
if (params.help) {
    log.info ''
    log.info 'G R A P E ~ RNA-seq Pipeline'
    log.info '----------------------------'
    log.info 'Run the GRAPE RNA-seq pipeline on a set of data.'
    log.info ''
    log.info 'Usage: '
    log.info '    grape-pipeline.nf -i INDEX_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...'
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --index INDEX_FILE                  Index file.'
    log.info '    --genome GENOME_FILE                Reference genome file(s).'
    log.info '    --annotation ANNOTAION_FILE         Reference gene annotation file(s).'
    log.info '    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all'
//    log.info '    --chunk-size CHUNK_SIZE             The number of records to be put in each chunk when splitting the input. Default: no split'
    log.info '    --error-strategy ERROR_STRATEGY     Specify how an error condition is managed by the pipeline processes. Possible values: ignore, retry'
    log.info '                                        Default: the entire pipeline  terminates if a process returns an error status.'
    log.info '    --max-mismatches THRESHOLD          Set maps with more than THRESHOLD error events to unmapped. Default "4".'
    log.info '    --max-multimaps THRESHOLD           Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10".'
    log.info ''
    log.info 'SAM read group options:'
    log.info '    --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.'
    log.info '    --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.'
    log.info '    --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.'
    log.info '    --rg-desc DESCRIPTION               Description for the BAM @RG tag.'
//    log.info '    --loglevel LOGLEVEL                 Log level (error, warn, info, debug). Default "info".'
    log.info ''
    exit 1
}

// check mandatory options
if (!params.genome) {
    exit 1, "Genome file not specified"
}

if (!params.annotation) {
    exit 1, "Annotation file not specified"
}

log.info ""
log.info "G R A P E ~ RNA-seq Pipeline"
log.info ""
log.info "General parameters"
log.info "------------------"
log.info "Index file                      : ${params.index}"
log.info "Genome                          : ${params.genome}"
log.info "Annotation                      : ${params.annotation}"
log.info "Pipeline steps                  : ${pipelineSteps.join(" ")}"
//log.info "Input chunk size                : ${params.chunkSize != null ? params.chunkSize : 'no split'}"
log.info "Error strategy                  : ${params.errorStrategy != null ? params.errorStrategy : 'default'}"
log.info ""

if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Max mismatches                  : ${params.maxMismatches}"
    log.info "Max multimaps                   : ${params.maxMultimaps}"
    log.info "Produce BAM stats               : ${params.bamStats}"
    if ( params.rgPlatform ) log.info "Sequencing platform             : ${params.rgPlatform}"  
    if ( params.rgLibrary ) log.info "Sequencing library              : ${params.rgLibrary}"  
    if ( params.rgCenterName ) log.info "Sequencing center               : ${params.rgCenterName}" 
    if ( params.rgDesc ) log.info "@RG Descritpiton                : ${params.rgDesc}" 
    log.info ""
}
if ('bigwig' in pipelineSteps) {
    log.info "Bigwig parameters"
    log.info "-----------------"
    log.info "References prefix               : ${params.wigRefPrefix != null ? params.wigRefPrefix : 'all'}"
    log.info ""
}

//if ('quantification' in pipelineSteps) {
//    log.info "Quantification parameters"
//    log.info "-------------------------"
//    log.info "Additional quantified elements  : ${params.countElements.size()==0 ? 'NONE' : params.countElements.join(" ")}"
//    log.info "Memory                          : ${params.fluxMem}"
//    log.info "Create profile file             : ${params.fluxProfile}"
//    log.info ""
//}

genomes=params.genome.split(',').collect { file(it) }
annos=params.annotation.split(',').collect { file(it) }

index = params.index ? file(params.index) : System.in
input_files = Channel.create()
input_chunks = Channel.create()

data = ['samples': [], 'ids': []]
merge = false
input = Channel
    .from(index.readLines())
    .map {
        line -> [ line.split()[0], line.split()[1], file(line.split()[2]), line.split()[3], line.split()[4] ]
    }

if (params.chunkSize) {
    merge = true
    input = input.splitFastq(by: params.chunkSize, file: true, elem: 2)
}

input.subscribe onNext: { 
        sample, id, path, type, view ->
        items = "sequencing runs"
        if( params.chunkSize ) {
            items = "chunks         "
            id = id+path.baseName.find(/\..+$/)
        } 
        input_chunks << tuple(sample, id, path, type, view)
        data['samples'] << sample
        data['ids'] << id }, 
    onComplete: { 
        ids=data['ids'].unique().size()
        samples=data['samples'].unique().size()
        if (ids != samples) merge=true
        log.info "Dataset information"
        log.info "-------------------"
        log.info "Number of sequenced samples     : ${samples}"
        log.info "Number of ${items}       : ${ids}" 
        log.info "Merging                         : ${ merge ? 'by sample' : 'none' }"
        log.info ""
        input_chunks << Channel.STOP 
    }

input_bams = Channel.create()
bam = Channel.create()

input_chunks
    .groupBy {
        sample, id, path, type, view -> id 
    }
    .flatMap ()
    .choice(input_files, input_bams) { it -> if ( it.value[0][3] == 'fastq' ) 0 else if ( it.value[0][3] == 'bam' ) 1}
    
input_files = input_files.map {        
    [it.key, it.value[0][0], it.value.collect { sample, id, path, type, view -> path }, fastq(it.value[0][2]).qualityScore()]
}

// id, sample, type, view, "${id}${prefix}.bam", pairedEnd
input_bams.map {
    [it.key, it.value[0][0], it.value[0][3], it.value[0][4], it.value.collect { sample, id, path, type, view -> path }, true].flatten()
}
.subscribe onNext: {
    bam << it
}, onComplete: {}

Genomes = Channel.create()
Annotations = Channel.create()
Channel.from(genomes)
    .merge(Channel.from(annos)) {
        g, a -> [g.name.split("\\.", 2)[0], g, a]
    }
    .separate(Genomes, Annotations) {
        a -> [[a[0],a[1]], [a[0],a[2]]]
    }

(Genomes1, Genomes2, Genomes3) = Genomes.into(3)
(Annotations1, Annotations2, Annotations3, Annotations4, Annotations5, Annotations6, Annotations7) = Annotations.into(7)

pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

if ('contig' in pipelineSteps || 'bigwig' in pipelineSteps) {
    process fastaIndex {

        input:
        set species, file(genome) from Genomes1
        set species, file(annotation) from Annotations1
    
        output:
        set species, file("${genome}.fai") into FaiIdx
    
        script:
        template(task.command)

    }
} else {
    FaiIdx = Channel.just(Channel.STOP)
}

(FaiIdx1, FaiIdx2) = FaiIdx.into(2)
 
if ('mapping' in pipelineSteps) {

    process index {
    
        input:
        set species, file(genome) from Genomes2
        set species, file(annotation) from Annotations7
    
        output:
        set species, file("genomeDir") into GenomeIdx
    
        script:
        sjOverHang = params.sjOverHang

        template(task.command)

    }
     
    (GenomeIdx1, GenomeIdx2) = GenomeIdx.into(2)

    process mapping {
    
        input:
        set id, sample, file(reads), qualityOffset from input_files
        set species, file(annotation) from Annotations3.first()
        set species, file(genomeDir) from GenomeIdx2.first()
    
        output:
        set id, sample, type, view, file("*.bam"), pairedEnd into bam
    
        script:
        type = 'bam'
        view = 'Alignments'
        prefix = pref
        maxMultimaps = params.maxMultimaps
        maxMismatches = params.maxMismatches
        
        // prepare BAM @RG tag information
        // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
        date = ""
        readGroup = []
        readGroup << "ID:${id}" 
        readGroup << "PU:${id}" 
        readGroup << "SM:${sample}" 
        if ( date ) readGroup << "DT:${date}"
        if ( params.rgPlatform ) readGroup << "PL:${params.rgPlatform}"
        if ( params.rgLibrary ) readGroup << "LB:${params.rgLibrary}"
        if ( params.rgCenterName ) readGroup << "CN:${params.rgCenterName}"
        if ( params.rgDesc ) readGroup << "DS:${params.rgDesc}"
        readGroup = readGroup.join(" ")

        fqs = reads.toString().split(" ")
        pairedEnd = false
        if (fqs.size() == 2) pairedEnd = true 

        template(task.command)
   
    }

    bam = bam.flatMap  { id, sample, type, view, path, pairedEnd ->
        [path].flatten().collect { f ->
            [id, sample, type, (f.name =~ /toTranscriptome/ ? 'Transcriptome' : 'Genome') + view, f, pairedEnd]
        }
    }

} else {
    GenomeIdx = Channel.just(Channel.STOP)
}

if ('quantification' in pipelineSteps) {
    
    process t_index {
    
        input:
        set species, file(genome) from Genomes3
        set species, file(annotation) from Annotations2

        output:
        set species, file('txDir') into TranscriptIdx
    
        script:
        template(task.command)

    }

} else {
    TranscriptIdx = Channel.just(Channel.STOP)
}

    
if (merge) {

    singleBam = Channel.create()
    groupedBam = Channel.create()
   
    bam.groupTuple(by: [1, 2, 3, 5]) // group by sample, type, view, pairedEnd (to get unique values for keys)
    .choice(singleBam, groupedBam) {
      it[4].size() > 1 ? 1 : 0
    }
    
    process mergeBam {
        
        input:
        set id, sample, type, view, file(bam), pairedEnd from groupedBam
    
        output:
        set id, sample, type, view, outfile, pairedEnd into mergedBam
    
        script:
        prefix = pref
        id = id.sort().join(':')
        outfile = "${sample}${prefix}_to${view.replace('Alignments','')}.bam"

        template(task.command)

    }
    
    bam = singleBam
    .mix(mergedBam)
    .map { 
        it.flatten() 
    }
}

if (!('mapping' in pipelineSteps)) {
    bam << Channel.STOP
}

(bam1, bam2) = bam
.into(2)

process inferExp {
    input:
    set id, sample, type, view, file(bam), pairedEnd from bam1
    set species, file(annotation) from Annotations4.first()

    output:
    set id, stdout into bamStrand

    script:
    prefix = pref
    genePred = "${annotation.name.split('\\.', 2)[0]}.genePred"
    bed12 = "${annotation.name.split('\\.', 2)[0]}.bed"

    template(task.command)
}


allBams = bam2.mix(bamStrand)
.groupBy()
.flatMap {
    bam -> bam.collect { [ it.value[0], it.value[1][-1] ].flatten() } 
}

(bam1, bam2, bam3, out) = allBams.into(4)

if (!('bigwig' in pipelineSteps)) bam1 = Channel.just(Channel.STOP)
if (!('contig' in pipelineSteps)) bam2 = Channel.just(Channel.STOP)
if (!('quantification' in pipelineSteps)) bam3 = Channel.just(Channel.STOP)

process bigwig {
    
    input:
    set id, sample, type, view, file(bam), pairedEnd, readStrand from bam1
    set species, file(genomeFai) from FaiIdx1.first()
    
    output:
    set id, type, views, file('*.bw') into bigwig

    script:
    type = "bigwig"
    wigRefPrefix = params.wigRefPrefix ?: ""
    views = task.views
    
    template(task.command)

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
    set id, sample, type, view, file(bam), pairedEnd, readStrand from bam2
    set species, file(genomeFai) from FaiIdx2.first()

    output:
    set id, type, view, file('*_contigs.bed') into contig

    script:
    type = 'bed'
    view = 'Contigs'

    template(task.command)

}

process quantification {

    input:    
    set id, sample, type, view, file(bam), file(bai), pairedEnd, readStrand from bam3.map { [it[0],it[1],it[2],file("${it[3].toAbsolutePath().toString().replace('.bam','.toTranscriptome.bam')}"),file("${it[3].toAbsolutePath()}.bai"),it[4],it[5]] }
    set species, file(txDir) from TranscriptIdx.first()

    output:
    set id, type, viewTx, file("Quant.genes.results") into isoforms
    set id, type, viewGn, file("Quant.isoforms.results") into genes

    script:
    type = "gtf"
    viewTx = "Transcript${txDir.name.replace('.gtf','').capitalize()}"
    viewGn = "Gene${txDir.name.replace('.gtf','').capitalize()}"
    memory = task.memory.toMega()
    
    template(task.command)
}

out.mix(bigwig, contig, isoforms, genes).collectFile(name: "pipeline.db", newLine: true) {
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
