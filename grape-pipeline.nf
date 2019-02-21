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

//Set default values for params
params.addXs = false
params.bamSort = null
params.chunkSize = null
params.dbFile = 'pipeline.db'
params.genomeIndex = null
params.help = false
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
params.wigRefPrefix = 'chr'

require 'modules.nf', params:params

// Some configuration variables
mappingTool = "${config.process.'withName:mapping'.ext.tool} ${config.process.'withName:mapping'.ext.version}"
bigwigTool = "${config.process.'withName:bigwig'.ext.tool} ${config.process.'withName:bigwig'.ext.version}"
quantificationTool = "${config.process.'withName:quantification'.ext.tool} ${config.process.'withName:quantification'.ext.version}"
quantificationMode = "${config.process.'withName:quantification'.ext.mode}"
useContainers = config.docker?.enabled ? 'docker' : (config.singularity?.enabled ? 'singularity' : 'no')
errorStrategy = config.process.errorStrategy
executor = config.process.executor ?: 'local'
queue = config.process.queue

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
//log.info "Input chunk size                : ${params?.chunkSize ?: 'no split'}"
log.info ""

if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Tool                            : ${mappingTool}"
    log.info "Max mismatches                  : ${params.maxMismatches}"
    log.info "Max multimaps                   : ${params.maxMultimaps}"
    //log.info "Produce BAM stats               : ${params.bamStats}"
    if ( params.rgPlatform ) log.info "Sequencing platform             : ${params.rgPlatform}"
    if ( params.rgLibrary ) log.info "Sequencing library              : ${params.rgLibrary}"
    if ( params.rgCenterName ) log.info "Sequencing center               : ${params.rgCenterName}"
    if ( params.rgDesc ) log.info "@RG Descritpiton                : ${params.rgDesc}"
    log.info ""
}
if ('bigwig' in pipelineSteps) {
    log.info "Bigwig parameters"
    log.info "-----------------"
    log.info "Tool                            : ${bigwigTool}"
    log.info "References prefix               : ${params?.wigRefPrefix ?: 'all'}"
    log.info ""
}

if ('quantification' in pipelineSteps) {
    log.info "Quantification parameters"
    log.info "-------------------------"
    log.info "Tool                            : ${quantificationTool}"
    log.info "Mode                            : ${quantificationMode}"
    log.info ""
}

log.info "Execution information"
log.info "---------------------"
log.info "Engine                          : ${executor}"
if (queue && executor != 'local')
    log.info "Queue(s)                        : ${queue}"
log.info "Use containers                  : ${useContainers}"
log.info "Error strategy                  : ${errorStrategy}"
log.info ""

genomes = params.genome.split(',').collect { file(it) }
annos = params.annotation.split(',').collect { file(it) }
if (params.genomeIndex) {
    genomeidxs = params.genomeIndex.split(',').collect { file(it) }
}

index = params.index ? file(params.index) : System.in
input_files = Channel.create()
input_chunks = Channel.create()

data = ['samples': [], 'ids': []]
input = Channel
    .from(index.readLines())
    .filter { it }  // get only lines not empty
    .map { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        return [sampleId, runId, resolveFile(fileName, index), format, readId]
    }

if (params.chunkSize) {
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
        log.info "Dataset information"
        log.info "-------------------"
        log.info "Number of sequenced samples     : ${samples}"
        log.info "Number of ${items}       : ${ids}"
        log.info "Merging                         : ${ ids != samples ? 'by sample' : 'none' }"
        log.info ""
        input_chunks << Channel.STOP
    }


workflow {

    flat_chunks = input_chunks.groupTuple(by:[0,1,3], sort:true)
    input_files = flat_chunks.filter { it[3] == 'fastq' }.map { sample, id, reads, type, view -> tuple( sample, id, reads, fastq(reads[0]).qualityScore() ) }
    input_bams = flat_chunks.filter { it[3] == 'bam' }

    input_bams.map {
        it.value.collect { bam ->
            [it.key, bam[0], bam[3], bam[4], bam[2], params.pairedEnd]
        }
    }.flatMap ()
    .subscribe onNext: {
        bam << it
    }, onComplete: {}

    def msg = "Output files db"
    log.info "=" * msg.size()
    log.info msg + " -> ${pdb}"
    log.info "=" * msg.size()
    log.info ""


    Genomes = Channel.create()
    Annotations = Channel.create()
    Channel.from(genomes)
        .merge(Channel.from(annos)) {
            g, a -> [g.name.split("\\.", 2)[0], g, a]
        }
        .separate(Genomes, Annotations) {
            a -> [[a[0],a[1]], [a[0],a[2]]]
        }

    FaiIdx = (  ('contig' in pipelineSteps || 'bigwig' in pipelineSteps) 
                ? fastaIndex( Genomes, Annotations ) 
                : Channel.empty() )

    if ('mapping' in pipelineSteps) {

        if (! params.genomeIndex) {
            GenomeIdx = index( Genomes, Annotations ) 

        } else {
            GenomeIdx = Genomes.merge(Channel.from(genomeidxs)) { d, i -> 
                [d[0], i]
            }
        }

        mapping( input_files, Annotations.first(), GenomeIdx.first() )

        bam = mapping.output.flatMap  { id, sample, type, view, path, pairedEnd ->
            [path].flatten().collect { f ->
                [id, sample, type, (f.name =~ /toTranscriptome/ ? 'Transcriptome' : 'Genome') + view, f, pairedEnd]
            }
        }

    } else {
        GenomeIdx = Channel.empty()
    }

    if ('quantification' in pipelineSteps && quantificationMode != "Genome") {

        QuantificationRef = txIndex(Genomes, Annotations)

    } else {
        QuantificationRef = Annotations
    }

    // Merge or rename bam

    groupped_bam = bam.groupTuple(by: [1, 2, 3, 5]) // group by sample, type, view, pairedEnd (to get unique values for keys)
    singleBam = groupped_bam.filter { it[4].size() <= 1 } 
    groupedBam = groupped_bam.filter { it[4].size() > 1 } 

    mergeBam( groupedBam )

    bam = singleBam
        .mix(mergeBam.output)
        .map {
            it.flatten()
        }

    if (!('mapping' in pipelineSteps)) {
        bam << Channel.STOP
    }

    if (params.readStrand) {
        
        bam.filter { it[3] =~ /Genome/ }
            .map {
                [it[0], params.readStrand]
            }.set { bamStrand }

    } else {

        bamStrand = inferExp(
            bam.filter { it[3] =~ /Genome/ },
            Annotations.first()
        )

    }

    allBams = bamStrand.cross(bam)
        .map {
            bam -> bam[1].flatten() + [bam[0][1]]
        }

    contigBams = bigwigBams = allBams.filter { it[3] =~ /Genome/ }
    quantificationBams = allBams.filter { it[3] =~ /${quantificationMode}/ }

    if (!('bigwig' in pipelineSteps)) bigwigBams = Channel.empty()
    if (!('contig' in pipelineSteps)) contigBams = Channel.empty()
    if (!('quantification' in pipelineSteps)) quantificationBams = Channel.empty()

    bigwig( bigwigBams, FaiIdx.first() ) 

    bigwig = bigwig.output.reduce([:]) { files, tuple ->
        def (id1, sample1, type1, view1, path1, pairedEnd1, readStrand1) = tuple
        if (!files) files = []
        paths = path1.toString().replaceAll(/[\[\],]/,"").split(" ").sort()
        (1..paths.size()).each { files << [id1, sample1, type1, view1[it-1], paths[it-1], pairedEnd1, readStrand1] }
        return files
    }
    .flatMap()

    contig( contigBams, FaiIdx.first())

    quantification( quantificationBams,  QuantificationRef.first() )

    allBams
        .mix(bigwig, contig.output, quantification.output.first, quantification.output.second )
        .collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { id, sample, type, view, file, pairedEnd, readStrand ->
            [sample, id, file, type, view, pairedEnd ? 'Paired-End' : 'Single-End', readStrand].join("\t")
        }
        .subscribe {
            log.info ""
            log.info "-----------------------"
            log.info "Pipeline run completed."
            log.info "-----------------------"
        }

}