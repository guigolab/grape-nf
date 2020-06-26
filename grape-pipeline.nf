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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

// imports
import groovy.json.JsonSlurper

//Set default values for params
params.addXs = false
params.mappingSortTool = null
params.chunkSize = null
params.dbFile = 'pipeline.db'
params.genomeIndex = null
params.help = false
params.markDuplicates = false
params.removeDuplicates = false
params.maxMismatches = 4
params.maxMultimaps = 10
params.pairedEnd = false
params.readLength = 150
params.readStrand = null
params.rgCenterName = null
params.rgDesc = null
params.rgLibrary = null
params.rgPlatform = null
params.sjOverHang = 100
params.steps = 'mapping,bigwig,contig,quantification'
params.wigRefPrefix = ''
params.inferExpThreshold = 0.8

// Process channels
Channel.empty().into {
    fetchInput
    fastaIndexGenomes; fastaIndexAnnotations
    indexRefs
    mappingInput; mappingIndex; mappingAnnotations
    txIndexRefs
    mergeBamInput
    markdupInput
    inferExpInput; inferExpAnnotations
    bamStatsInput; bamStatsAnnotations
    bigwigInput; bigwigFastaIndex
    contigInput; contigFastaIndex
    quantificationInput; quantificationTxIndex
}

// Auxiliary variables
def comprExts = ['gz', 'bz2', 'zip']
def pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

// Clear pipeline.db file
pdb = file(params.dbFile)
pdb.write('')

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
    log.info '    grape-pipeline.nf --index INDEX_FILE --genome GENOME_FILE --annotation ANNOTATION_FILE [OPTION]...'
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --index INDEX_FILE                  Index file.'
    log.info '    --genome GENOME_FILE                Reference genome file(s).'
    log.info '    --annotation ANNOTAION_FILE         Reference gene annotation file(s).'
    log.info '    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all'
//    log.info '    --chunk-size CHUNK_SIZE             The number of records to be put in each chunk when splitting the input. Default: no split'
    log.info '    --max-mismatches THRESHOLD          Set maps with more than THRESHOLD error events to unmapped. Default "4".'
    log.info '    --max-multimaps THRESHOLD           Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10".'
    log.info '    --bam-sort METHOD                   Specify the method used for sorting the genome BAM file.'
    log.info '    --paired-end                        Treat input data as paired-end.'
    log.info '    --add-xs                            Add the XS field required by Cufflinks/Stringtie to the genome BAM file.'
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
if (!params.genomeIndex && !params.genome) {
    exit 1, "Reference genome not specified"
}

if ('quantification' in pipelineSteps && !params.annotation) {
    exit 1, "Annotation not specified"
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
log.info "Pipeline profile                : ${workflow.profile}"
log.info ""

if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    // log.info "Tool                            : ${mappingTool}"
    log.info "Max mismatches                  : ${params.maxMismatches}"
    log.info "Max multimaps                   : ${params.maxMultimaps}"
    if ( params.rgPlatform ) log.info "Sequencing platform             : ${params.rgPlatform}"
    if ( params.rgLibrary ) log.info "Sequencing library              : ${params.rgLibrary}"
    if ( params.rgCenterName ) log.info "Sequencing center               : ${params.rgCenterName}"
    if ( params.rgDesc ) log.info "@RG Descritpiton                : ${params.rgDesc}"
    log.info ""
}
if ('bigwig' in pipelineSteps) {
    log.info "Bigwig parameters"
    log.info "-----------------"
    // log.info "Tool                            : ${bigwigTool}"
    log.info "References prefix               : ${['','-'].contains(params.wigRefPrefix) ? 'all' : params.wigRefPrefix}"
    log.info ""
}

if ('quantification' in pipelineSteps) {
    log.info "Quantification parameters"
    log.info "-------------------------"
    log.info "Tool                            : ${params.quantificationTool}"
    log.info "Mode                            : ${params.quantificationMode}"
    log.info ""
}

log.info "Execution information"
log.info "---------------------"
log.info "Use containers                  : ${workflow.containerEngine?:false}"
log.info ""

def msg = "Output files db"
log.info "=" * msg.size()
log.info msg + " -> ${pdb}"
log.info "=" * msg.size()
log.info ""

// Get input data
index = params.index ? file(params.index) : System.in

Channel.from(index.readLines())
.filter { it }  // get only non-empty lines
.map { line ->
    def (sampleId, runId, fileName, format, readId) = line.split()
    def fetch = false
    if ( fileName.split(',').size() > 1 )
        fetch = true
    if ( ! fetch )
        fileName = resolveFile(fileName, index)
    [sampleId, runId, fileName, format, readId, fetch]
}.tap {
    inputFilesForFetch
    inputFiles
}

(ids, samples, indexLines) = readTsv(index)

log.info "Dataset information"
log.info "-------------------"
log.info "Number of sequenced samples     : ${samples}"
log.info "Number of sequencing runs       : ${ids}"
log.info "Merging                         : ${ ids != samples ? 'by sample' : 'none' }"
log.info ""

inputFilesForFetch
.filter { it[5] }
.map { sampleId, runId, fileName, format, readId, fetch ->
    [sampleId, runId, fileName, format, readId]
}.set {
    fetchInput
}

inputFiles.filter { !it[5] }
.map { sampleId, runId, fileName, format, readId, fetch ->
    [sampleId, runId, fileName, format, readId]
}.set {
    inputFilesNotToFetch
}

// // Get references
// genomes = file(params.genome)
// genomes.sort()
// if (params.oneGenome) {
//    genomes = [genomes]
// }
// annotations = file(params.annotation)
// if (params.genomeIndex) {
//     genomeidxs = file(params.genomeIndex)
// }

// Channel.from(genomes)
// .combine(Channel.from(annotations)).view()
// .into { refsForGenomes; refsForAnnotations }

// refsForGenomes.filter {
//     it[1].name =~ /.fa(.gz)?$/
// }.set{ Genomes }

// refsForAnnotations.filter {
//     it[1].name =~ /.gtf(.gz)?$/
// }.set{ Annotations }

Channel.fromPath(params.genome).toSortedList().flatten().into {
    genomesForFastaIndex
    genomesForIndex
    genomesForTxIndex
}
Channel.fromPath(params.annotation).toSortedList().flatten().into {
    annotationsForFastaIndex
    annotationsForIndex
    annotationsForMapping
    annotationsForTxIndex
    annotationsForInferExp
    annotationsForBamStats
    annotationsForQuantification
}

if ( 'bigwig' in pipelineSteps || 'contig' in pipelineSteps) {
    genomesForFastaIndex.set { fastaIndexGenomes }
    annotationsForFastaIndex.set { fastaIndexAnnotations }
}

if ( 'mapping' in pipelineSteps ) {
    if ( params.genomeIndex ) {
        Channel.from(genomeidxs).map {
            [ it.baseName, it ]
        }.set { mappingIndex }
    } else {
        // genomesForIndex.set { indexGenomes }
        // annotationsForIndex.set { indexAnnotations }
        genomesForIndex.merge(annotationsForIndex).set { indexRefs }
    }
    annotationsForMapping.set { mappingAnnotations }
    annotationsForInferExp.set { inferExpAnnotations }
    annotationsForBamStats.set { bamStatsAnnotations }
}

if ( 'quantification' in pipelineSteps && params.quantificationMode == "Transcriptome" ) {
    // genomesForTxIndex.set { txIndexGenomes }
    // annotationsForTxIndex.set { txIndexAnnotations }
    genomesForTxIndex.merge(annotationsForTxIndex).set { txIndexRefs }
}

// Processes
process fetch {
    tag { outPath.name }
    storeDir { outPath.parent }

    input:
    set sample, id, path, type, view from fetchInput

    output:
    set sample, id, file("${outPath.name}"), type, view into fetchOutput

    script:
    def paths = path.split(',')
    outPath = workflow.launchDir.resolve(paths[-1])
    urls = paths.size() > 1 ? paths[0..-2].join(' ') : ''
    """
    for url in $urls; do
        if wget \${url}; then
            exit 0
        fi
    done
    """
}

inputFilesNotToFetch.mix(fetchOutput)
.groupTuple(by: [0,1,3], sort: true)
.into {
    inputFilesForBams
    inputFilesForFastqs
}

inputFilesForFastqs.filter {
    it[3] == 'fastq'
}.map {
    [it[1], it[0], it[2], fastq(it[2][0]).qualityScore()]
}.set { mappingInput }

inputFilesForBams.filter {
    it[3] == 'bam'
}.transpose()
.map { sample, id, path, type, view ->
    def genome = file(params.genome)
    def annotation = file(params.annotation)
    if ( genome instanceof LinkedList || annotation instanceof LinkedList ) {
        log.error "Cannot start from bam with multiple references"
        exit(1)
    }
    def ref = "${genome.baseName}_${annotation.baseName}"
    [ref, id, sample, type, view, path, params.pairedEnd].flatten()
}
.set {
    inputBams
}

process fastaIndex {

    tag "${genome.baseName}-${params.fastaIndexTool}-${params.fastaIndexToolVersion}"

    input:
    file(genome) from fastaIndexGenomes

    output:
    file faiIndex into fastaIndexOutput

    script:
    def gz = genome.name.lastIndexOf('.gz') > -1 ?: genome.name.length()
    faiIndex = "${genome.name.substring(0, gz)}.fai"
    compressed = genome.extension in comprExts ? "-${genome.extension}" : ''
    command = "${task.process}/${params.fastaIndexTool}${compressed}"
    template(command)

}

fastaIndexOutput.into {
    fastaIndexForBigwig
    fastaIndexForContig
}

if ( 'bigwig' in pipelineSteps ) {
    fastaIndexForBigwig.set { bigwigFastaIndex }
}

if ( 'contig' in pipelineSteps ) {
    fastaIndexForContig.set { contigFastaIndex }
}

process index {

    label "mapping"
    tag "${prefix}-${params.mappingTool}-${params.mappingToolVersion}"

    input:
    // file(genome) from indexGenomes
    // file(annotation) from indexAnnotations
    set file(genome), file(annotation) from indexRefs

    output:
    set file(prefix), file(annotation) into indexOutput

    script:
    prefix = "${genome.baseName}_${annotation.baseName}"
    cpus = task.cpus
    sjOverHang = params.sjOverHang
    readLength = params.readLength
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${task.process}/${params.mappingTool}${genomeCompressed}${annoCompressed}"

    template(command)

}

if ( ! params.genomeIndex ) {
    indexOutput.set { mappingIndex }
}

process txIndex {

    label 'quantification'
    tag "${prefix}-${params.quantificationTool}-${params.quantificationToolVersion}"

    input:
    // file(genome) from txIndexGenomes
    // file(annotation) from txIndexAnnotations
    set file(genome), file(annotation) from txIndexRefs

    output:
    file(prefix) into txIndexOutput

    script:
    prefix = "${genome.baseName}_${annotation.baseName}"
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${task.process}/${params.quantificationTool}${genomeCompressed}${annoCompressed}"

    template(command)

}

process mapping {

    label "mapping"
    tag "${id.replace(':', '_')}-${params.mappingTool}-${params.mappingToolVersion}"

    input:
    set id, sample, file(reads), qualityOffset from mappingInput
    // each file(annotation) from mappingAnnotations
    each file(refs) from mappingIndex

    output:
    set ref, id, sample, type, view, file("*.bam"), pairedEnd into mappingOutput

    script:
    type = 'bam'
    view = 'Alignments'
    prefix = "${sample}${pref}"
    maxMultimaps = params.maxMultimaps
    maxMismatches = params.maxMismatches
    genomeDir = refs[0]
    annotation = refs[1]
    ref = genomeDir.name

    // prepare BAM @RG tag information
    // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
    date = ""
    readGroupList = []
    readGroupList << ["ID", "${id}"]
    readGroupList << ["PU", "${id}"]
    readGroupList << ["SM", "${sample}"]
    if ( date ) readGroupList << ["DT", "${date}"]
    if ( params.rgPlatform ) readGroupList << ["PL", "${params.rgPlatform}"]
    if ( params.rgLibrary ) readGroupList << ["LB", "${params.rgLibrary}"]
    if ( params.rgCenterName ) readGroupList << ["CN", "${params.rgCenterName}"]
    if ( params.rgDesc ) readGroupList << ["DS", "${params.rgDesc}"]
    (s,t) = params.mappingReadGroupSeparators
    readGroup = readGroupList.collect { it.join(s) }.join(t)

    fqs = reads.toString().split(" ")
    pairedEnd = (fqs.size() == 2)
    taskMemory = task.memory ?: 1.GB
    totalMemory = (taskMemory.toBytes()*2/3) as int
    threadMemory = (totalMemory/task.cpus) as int
    cpus = task.cpus
    halfCpus = (task.cpus > 1 ? task.cpus / 2 : task.cpus) as int

    command = "${task.process}/${params.mappingTool}-${params.mappingToolVersion.split("\\.")[0..1].join(".")}"
    switch(params.mappingTool) {
        case 'GEM':
            command += "-${pairedEnd ? 'Paired-End' : 'Single-End'}"
            break
        case 'STAR':
            command += (params.mappingSortTool ? "-"+params.mappingSortTool : '') + (params.quantificationMode ? "-"+params.quantificationMode : '') + (params.addXs ? "-XS" : '')
            break
    }
    template(command)

}

mappingOutput.flatMap  { ref, id, sample, type, view, path, pairedEnd ->
    [path].flatten().collect { f ->
        [ref, id, sample, type, (f.name =~ /toTranscriptome/ ? 'Transcriptome' : 'Genome') + view, f, pairedEnd]
    }
}.mix(inputBams).groupTuple(by: [0, 2, 3, 4, 6]) // group by sample, type, view, pairedEnd (to get unique values for keys)
.into {
    bamFilesForSingle
    bamFilesForMergeBam
}

bamFilesForSingle.filter {
    it[5].size() == 1
}.set {
    bamFilesSingle
}

bamFilesForMergeBam.filter {
    it[5].size() > 1
}.into {
    bamFilesGenomeForMerge
    bamFilesTranscriptomeForMerge
}

bamFilesGenomeForMerge.filter {
    it[4] =~ /^Genome/
}.set {
    mergeBamGenomeInput
}

bamFilesTranscriptomeForMerge.filter {
    it[4] =~ /^Transcriptome/
}.transpose()
.set {
    bamFilesTranscriptomeMerge
}

process sortBam {
    tag "${id}-${params.mergeBamTool}-${params.mergeBamToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd from bamFilesTranscriptomeMerge

    output:
    set ref, id, sample, type, view, file("${prefix}.bam"), pairedEnd into mergeBamTranscriptomeInput

    script:
    cpus = task.cpus
    taskMemory = task.memory ?: 1.GB
    totalMemory = taskMemory.toBytes()
    threadMemory = totalMemory/cpus
    prefix = "${bam.baseName}_sorted"
    command = "${task.process}/${params.mergeBamTool}"

    template(command)
}

mergeBamGenomeInput.mix(
    mergeBamTranscriptomeInput.groupTuple(by: [0, 2, 3, 4, 6], sort: true)
).set {
    mergeBamInput
}

process mergeBam {

    tag "${id.replace(':', '_')}-${params.mergeBamTool}-${params.mergeBamToolVersion}"

    input:
    set ref, id, sample, type, view, file("${sample}_??.bam"), pairedEnd from mergeBamInput

    output:
    set ref, id, sample, type, view, file("${prefix}.bam"), pairedEnd into mergeBamOutput

    script:
    cpus = task.cpus
    id = id.sort().join(':')
    prefix = "${sample}${pref}_to${view.replace('Alignments','')}"
    command = "${task.process}/${params.mergeBamTool}"

    template(command)

}

bamFilesSingle
.mix(mergeBamOutput)
.map {
    it.flatten()
}.into {
    bamFilesForGenome
    bamFilesForTranscriptome
}

bamFilesForTranscriptome.filter {
    it[4] =~ /^Transcriptome/
}.set{
    bamFilesToTranscriptome
}

bamFilesForGenome.filter {
    it[4] =~ /^Genome/
}.into{
    bamFilesForMarkdup
    bamFilesToGenome
}

if ( params.markDuplicates || params.removeDuplicates ) {
    bamFilesForMarkdup.set { markdupInput }
}

process markdup {

    tag "${id.replace(':', '_')}-${params.markdupTool}-${params.markdupToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd from markdupInput

    output:
    set ref, id, sample, type, view, file("${prefix}.bam"), pairedEnd into markdupOutput

    script:
    cpus = task.cpus
    memory = (task.memory ?: 2.GB).toMega()
    prefix = "${bam.baseName}.markdup"

    command = "${task.process}/${params.markdupTool}${params.removeDuplicates ? '-remove' : ''}"
    template(command)
}

if ( params.markDuplicates || params.removeDuplicates ) {
    bamFilesToGenome = markdupOutput
}

bamFilesToGenome.map {
    [ it[0].tokenize('_')[1] ] + it
}.combine(
    inferExpAnnotations.map {
        [it.baseName, it]
    }, by: 0
).map {
    it[1..-2] + [it[-1]]
}
.into {
    inferExpInput
    bamStatsInput
}

process inferExp {

    tag "${id.replace(':', '_')}-${params.inferExpTool}-${params.inferExpToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd, file(annotation) from inferExpInput

    output:
    // set id, stdout into bamStrand
    set ref, id, sample, type, view, file(bam), pairedEnd, stdout into inferExpOutputJSON

    script:
    prefix = "${annotation.name.split('\\.', 2)[0]}"
    command = "${task.process}/${params.inferExpTool}"
    threshold = params.inferExpThreshold

    template(command)
}

inferExpOutputJSON.map {
    j = new JsonSlurper()
    d = j.parseText(it[-1])
    it[0..-3] + [ d.paired, d.exp ]
}.set {
    inferExpOutput
}

inferExpOutput.into {
    bigwigContigInput
    bamFilesCrossTranscriptome
    bamFilesCrossBamStats
    quantificationInputGenome
    bamFilesToGenome
}

bamFilesToTranscriptome.combine(bamFilesCrossTranscriptome, by: [0,1]).map {
    it[0..5] + it[-2..-1]
}.into {
    bamFilesToTranscriptome
    quantificationInputTranscriptome
}

switch(params.quantificationMode) {
    case 'Genome':
        quantificationInputGenome.map {
            [ it[0].tokenize('_')[1] ] + it
        }.combine(
            annotationsForQuantification.map {
                [it.baseName, it]
            }, by: 0
        ).map {
            it[1..-2] + [it[-1]]
        }
        .set {
            quantificationInput
        }
        break
    case 'Transcriptome':
        quantificationInputTranscriptome.combine(
            txIndexOutput.map {
                [it.name, it]
            }, by: 0
        )
        .map {
            it[0..-2] + [it[-1]]
        }
        .set {
            quantificationInput
        }
        break
}

process bamStats {

    tag "${id.replace(':', '_')}-${params.bamStatsTool}-${params.bamStatsToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd, file(annotation) from bamStatsInput

    output:
    set ref, id, sample, type, views, file('*.json'), pairedEnd into bamStatsOutput

    script:
    cpus = task.cpus
    type = "json"
    prefix = "${sample}"
    views = "BamStats"
    maxBuf = params.bamStatsMaxBuf
    logLevel = params.bamStatsLogLevel
    command = "${task.process}/${params.bamStatsTool}"

    template(command)
}

bamStatsOutput.combine(bamFilesCrossBamStats, by: [0,1]).map {
    it[0..6] + it[-1]
}.set {
    bamStatsFiles
}

bigwigContigInput.map {
    [ it[0].tokenize('_')[0] ] + it
}.combine(
    bigwigFastaIndex.map {
        def name = it.baseName
        [name.substring(0, name.lastIndexOf('.fa')) , it]
    }, by: 0
).map {
    it[1..-2] + [it[-1]]
}
.into {
    bigwigInput
    contigInput
}

process bigwig {

    tag "${id.replace(':', '_')}-${params.bigwigTool}-${params.bigwigToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd, readStrand, file(genomeFai) from bigwigInput

    output:
    set ref, id, sample, type, views, file('*.bw'), pairedEnd, readStrand into bigwigOutput

    script:
    cpus = task.cpus
    type = "bigWig"
    prefix = "${sample}"
    wigRefPrefix = params.wigRefPrefix != "-" ? params.wigRefPrefix : ""
    views = params.bigwigViews[readStrand]
    command = "${task.process}/${params.bigwigTool}-${readStrand}"

    template(command)

}

process contig {

    tag "${id.replace(':', '_')}-${params.contigTool}-${params.contigToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd, readStrand, file(genomeFai) from contigInput

    output:
    set ref, id, sample, type, view, file('*.bed'), pairedEnd, readStrand into contigOutput

    script:
    cpus = task.cpus
    type = 'bed'
    view = 'Contigs'
    prefix = "${sample}.contigs"
    command = "${task.process}/${params.contigTool}-${readStrand}"

    template(command)

}

process quantification {

    label 'quantification'
    tag "${id.replace(':', '_')}-${params.quantificationTool}-${params.quantificationToolVersion}"

    input:
    set ref, id, sample, type, view, file(bam), pairedEnd, readStrand, file(quantRef) from quantificationInput

    output:
    set ref, id, sample, type, viewTx, file("*isoforms*"), pairedEnd, readStrand into quantificationIsoforms
    set ref, id, sample, type, viewGn, file("*genes*"), pairedEnd, readStrand into quantificationGenes

    script:
    cpus = task.cpus
    prefix = "${sample}"
    type = params.quantificationFileType
    viewTx = "Transcript"
    viewGn = "Gene"
    memory = (task.memory ?: 1.GB).toMega()
    command = "${task.process}/${params.quantificationTool}"
    if ( params.quantificationTool == 'RSEM') {
        command += "-${pairedEnd ? 'Paired-End' : 'Single-End'}"
    }
    command += "-${readStrand}"

    template(command)

}

bigwigOutput.flatMap { ref, id, sample, type, views, files, pairedEnd, readStrand ->
    [views, files].transpose().collect { view, f ->
        [ ref, id, sample, type, view, f, pairedEnd, readStrand ]
    }
}.set {
    bigwigFiles
}

bamFilesToGenome.mix(bamFilesToTranscriptome, bamStatsFiles, bigwigFiles, contigOutput, quantificationIsoforms, quantificationGenes)
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { ref, id, sample, type, view, file, pairedEnd, readStrand ->
    if ( !view.contains(ref) ) {
        r = ref.tokenize('_')
        view = "${view}${r.collect { it.capitalize() }.join('')}"
    }
    [ sample, id, file, type, view, pairedEnd ? 'Paired-End' : 'Single-End', readStrand].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}

/*
 * Given the input index file returns the number of unique samples,
 * the number of unique runs, and the lines of the index.
 * Params:
 * - tsvFile: a file object representing the TSV file
 */
def readTsv(tsvFile) {
    def (samples, ids, lines) = [[], [], []]
    tsvFile.eachLine { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        samples << sampleId
        ids << runId
        lines << line
    }
    [ids.unique().size(), samples.unique().size(), lines]
}

/*
 * Given a string path resolve it against the index file location.
 * Params:
 * - str: a string value represting the file pah to be resolved
 * - index: path location against which relative paths need to be resolved
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str)
  }
}

def testResolveFile() {
  def index = file('/path/to/index')
  assert resolveFile('str', index) == file('/path/to/str')
  assert resolveFile('/abs/file', index) == file('/abs/file')
  assert resolveFile('s3://abs/file', index) == file('s3://abs/file')
}
