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

nextflow.enable.dsl=2

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

Channel.empty().set { markdupInput }

// Import modules

include { fetch } from './modules/fetch.nf'
include { fastaIndex } from './modules/fastaindex.nf'
include { index } from './modules/index.nf'
include { txIndex } from './modules/txIndex.nf'
include { mapping } from './modules/mapping.nf'
include { sortBam } from './modules/sortBam.nf'
include { mergeBam } from './modules/mergeBam.nf'
include { inferExp } from './modules/inferExp.nf'
include { markdup } from './modules/markdup.nf'
include { contig } from './modules/contig.nf'
include { quantification } from './modules/quantification.nf'
include { bamStats } from './modules/bamStats.nf'
include { bigwig} from './modules/bigwig.nf'


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

// Get references
genomes = params.genome.split(',').collect { file(it) }
annotations = params.annotation.split(',').collect { file(it) }
if (params.genomeIndex) {
    genomeidxs = params.genomeIndex.split(',').collect { file(it) }
}

Channel.from(genomes)
.merge(Channel.from(annotations)) { g,a ->
    [ g.simpleName, [g, a] ]
}.transpose()
.set { refs }

refs.filter {
    it[1].name =~ /.fa(.gz)?$/
}.set{ Genomes }

refs.filter {
    it[1].name =~ /.gtf(.gz)?$/
}.set{ Annotations }



if ( 'bigwig' in pipelineSteps || 'contig' in pipelineSteps) {
    Genomes.set { fastaIndexGenomes }
    Annotations.set { fastaIndexAnnotations }
}

if ( 'mapping' in pipelineSteps ) {
    if ( params.genomeIndex ) {
        Channel.from(genomeidxs).map {
            [ it.simpleName, it ]
        }.set { mappingIndex
 }
    } else {
        Genomes.set { indexGenomes }
        Annotations.set { indexAnnotations }
    }
}

if ( 'quantification' in pipelineSteps && params.quantificationMode == "Transcriptome" ) {
    Genomes.set { txIndexGenomes }
    Annotations.set { txIndexAnnotations }
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


workflow {


  fetchOutput = fetch( fetchInput )

  inputFilesNotToFetch.mix(fetchOutput)
  .groupTuple(by: [0,1,3], sort: true)
  .set {
      inputFiles
  }


  inputFiles.filter {
      it[3] == 'fastq'
  }.map {
      [it[1], it[0], it[2], fastq(it[2][0]).qualityScore()]
  }.set { mappingInput }

  inputFiles.filter {
      it[3] == 'bam'
  }.transpose()
  .map { sample, id, path, type, view ->
      [id, sample, type, view, path, params.pairedEnd].flatten()
  }
  .set {
      inputBams
  }

  fastaIndexOutput = fastaIndex( fastaIndexGenomes , fastaIndexAnnotations )




  if ( 'bigwig' in pipelineSteps ) {
      fastaIndexOutput.set { bigwigFastaIndex }
  }

  if ( 'contig' in pipelineSteps ) {
      fastaIndexOutput.set { contigFastaIndex }
  }

  txIndexOutput = txIndex( txIndexGenomes , txIndexAnnotations )



  if ( ! params.genomeIndex ) {
      indexOutput = index( indexGenomes , indexAnnotations )
      indexOutput.set { mappingIndex }
  }

  mappingOutput = mapping( mappingInput , Annotations.first() , mappingIndex.first() )
  
  

  mappingOutput.flatMap  { id, sample, type, view, path, pairedEnd ->
    [path].flatten().collect { f ->
        [id, sample, type, (f.name =~ /toTranscriptome/ ? 'Transcriptome' : 'Genome') + view, f, pairedEnd]
    }
  }.mix(inputBams).groupTuple(by: [1, 2, 3, 5]) // group by sample, type, view, pairedEnd (to get unique values for keys)
  .set{
      bamFiles
  }
 
  

  bamFiles.filter {
      it[4].size() == 1
  }.set {
      bamFilesSingle
  }
  

  bamFiles.filter {
      it[4].size() > 1
  }.set {
       bamFilesForMerge
  }

  bamFilesForMerge.filter {
      it[3] =~ /^Genome/
  }.set {
      mergeBamGenomeInput
  }

  bamFilesForMerge.filter {
      it[3] =~ /^Transcriptome/
  }.transpose()
  .set {
      bamFilesTranscriptomeMerge
  }
  
  mergeBamTranscriptomeInput = sortBam( bamFilesTranscriptomeMerge )

   mergeBamGenomeInput.mix(
     mergeBamTranscriptomeInput.groupTuple(by: [1, 2, 3, 5], sort: true)
   ).set {
       mergeBamInput
   }

   mergeBamOutput = mergeBam( mergeBamInput )

  bamFilesSingle
  .mix(mergeBamOutput)
  .map {
      it.flatten()
  }.set{
      bamFilesForGenome
  }

  bamFilesForGenome.filter {
      it[3] =~ /^Genome/
  }.set{
      bamFilesForMarkdup
  }

  bamFilesForGenome.filter {
      it[3] =~ /^Genome/
  }.set{
      bamFilesToGenome
  }


  bamFilesSingle
  .mix(mergeBamOutput)
  .map {
      it.flatten()
  }.set{
      bamFilesForTranscriptome
  }

  bamFilesForTranscriptome.filter {
      it[3] =~ /^Transcriptome/
  }.set{
      bamFilesToTranscriptome
  }
  



  
  if ( params.markDuplicates || params.removeDuplicates ) {
     bamFilesForMarkdup.set { markdupInput }
  }

  
  markdupOutput = markdup( markdupInput )

  if ( params.markDuplicates || params.removeDuplicates ) {
    bamFilesToGenome = markdupOutput
  }
  
  inferExpOutputJSON = inferExp( bamFilesToGenome , Annotations.first() )

  
  
  inferExpOutputJSON.map {
    j = new JsonSlurper()
    d = j.parseText(it[-1])
    it[0..-3] + [ d.paired, d.exp ]
  }.set {
      inferExpOutput
  }



  bamFilesToTranscriptome.cross(inferExpOutput).map { transcriptome, genome ->
      transcriptome[0..-2] + genome[-2..-1]
  }.set {
      bamFilesToTranscriptome
  }

  switch(params.quantificationMode) {
      case 'Genome':
          inferExpOutput.set { quantificationInput }
          annotationsForQuantification.set { quantificationIndex }
          break
      case 'Transcriptome':
          bamFilesToTranscriptome.set { quantificationInput }
          txIndexOutput.set { quantificationIndex }
          break
  }


   bamStatsOutput = bamStats( bamFilesToGenome , Annotations.first() )

   bamStatsOutput.cross(inferExpOutput).map { stats, genome ->
     stats + genome[-1]
   }.set {
     bamStatsFiles
   }

   bigwigOutput = bigwig( inferExpOutput , bigwigFastaIndex.first() )
   

   contigOutput = contig( inferExpOutput , contigFastaIndex.first() )




   (quantificationIsoforms, quantificationGenes) = quantification( quantificationInput , quantificationIndex.first() )

   

  bigwigOutput.flatMap { id, sample, type, views, files, pairedEnd, readStrand ->
    [views, files].transpose().collect { view, f ->
        [ id, sample, type, view, f, pairedEnd, readStrand ]
    }
  }.set {
      bigwigFiles
  }




  inferExpOutput.mix(bamFilesToTranscriptome, bamStatsFiles, bigwigFiles, contigOutput, quantificationIsoforms, quantificationGenes)
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