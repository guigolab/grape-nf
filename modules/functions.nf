// imports
import groovy.json.JsonSlurper

// Print pipeline usage and exit if `--help` is passed
def printUsage() {
    if (params.help || params.helpAll) {
        log.info ''
        log.info 'G R A P E ~ RNA-seq Pipeline'
        log.info '----------------------------'
        log.info 'Run the GRAPE RNA-seq pipeline on a set of data.'
        log.info ''
        log.info 'Usage: '
        log.info '    grape-pipeline --index INDEX_FILE --genome GENOME_FILE --annotation ANNOTATION_FILE [OPTION]...'
        log.info ''
        log.info 'Options:'
        log.info '    --help                              Show this message and exit'
        log.info '    --index INDEX_FILE                  Index file'
        log.info '    --genome GENOME_FILE                Reference genome file'
        log.info '    --annotation ANNOTATION_FILE        Reference gene annotation file'
        log.info '    --genomeIndex GENOME_INDEX          Pre-existent genome index file for mapping (skip mapping index process)'
        log.info '    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all'
        log.info '    --maxMismatches THRESHOLD           Set maps with more than THRESHOLD error events to unmapped. Default "4"'
        log.info '    --maxMultimaps THRESHOLD            Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10"'
        log.info '    --readStrand DIRECTION              Specify the orientation of the reads for all samples'
        log.info '    --pairedEnd                         Treat input data as paired-end'
        log.info '    --addXs                             Add the XS field required by Cufflinks/Stringtie to the genome BAM file'
        log.info '    --sjOverHang                        The overhang length for junctions'
        log.info '    --markDuplicates                    Mark duplicates in the genomic alignemnt files'
        log.info '    --removeDuplicates                  Mark and remove duplicates from the genomic alignemnt files'
        log.info '    --inferExpThreshold                 Read percentage threshold for inferring library configuration'
        log.info ''
        log.info 'SAM read group options:'
        log.info '    --rgPlatform PLATFORM               Platform/technology used to produce the reads for the BAM @RG tag'
        log.info '    --rgLibrary LIBRARY                 Sequencing library name for the BAM @RG tag'
        log.info '    --rgCenterName CENTER_NAME          Name of sequencing center that produced the reads for the BAM @RG tag'
        log.info '    --rgDesc DESCRIPTION                Description for the BAM @RG tag'
        log.info ''
        log.info 'Output options:'
        log.info '    --dbFile DB_FILE                    Location of the pipeline DB file'
        if (params.helpAll) {
            log.info ''
            log.info 'Tool selection:'
            log.info '    --mappingTool                       Specify the mapping tool'
            // log.info '    --mappingToolVersion                Specify the mapping tool version'
            log.info '    --mappingSortTool                   Specify the sorting tool for the mapping step'
            log.info '    --bamStatsTool                      Specify the bamstats tool'
            log.info '    --bigwigTool                        Specify the bigwig tool'
            log.info '    --contigTool                        Specify the contig tool'
            log.info '    --fastaIndexTool                    Specify the fastaIndex tool'
            log.info '    --inferExpTool                      Specify the inferExp tool'
            log.info '    --markdupTool                       Specify the markdup tool'
            log.info '    --mergeBamTool                      Specify the mergeBam tool'
            log.info '    --quantificationTool                Specify the quantification tool'
            log.info ''
            log.info 'NOTE: To specify the version of a tool please use a camel case parameter starting with the tool name in lower case and then the word `version` (e.g. `--starVersion`)'
        }
        exit 1
    }
}

// Print pipeline log
def printLog() {
    def paramNames = params.keySet()
    
    def duplicatesMode = 'keep'
    if ( params.markDuplicates ) duplicatesMode = 'mark'
    if ( params.removeDuplicates ) duplicatesMode = 'remove'

    log.info ""
    log.info "G R A P E ~ RNA-seq Pipeline"
    log.info ""
    log.info "General parameters"
    log.info "------------------"
    log.info "Index file                      : ${file(params.index)}"
    log.info "Genome                          : ${file(params.genome)}"
    log.info "Annotation                      : ${file(params.annotation)}"
    log.info "Pipeline steps                  : ${params.stepList.join(" ")}"
    log.info "Pipeline profile                : ${workflow.profile}"
    log.info "Read strand                     : ${params.readStrand?:'unspecified'}"
    log.info "Paired-end                      : ${params.pairedEnd}"
    log.info "PCR duplicates mode             : ${duplicatesMode}"
    log.info "Use containers                  : ${workflow.containerEngine?:false}"
    log.info ""

    if ('mapping' in params.stepList) {
        log.info "Mapping parameters"
        log.info "------------------"
        log.info "Tool                            : ${params.mappingTool}"
        log.info "Max mismatches                  : ${params.maxMismatches}"
        log.info "Max multimaps                   : ${params.maxMultimaps}"
        log.info "Splice junctions overhang       : ${params.sjOverHang}"
        log.info "Add XS tag                      : ${params.addXs}"
        if ( params.rgPlatform ) log.info "Sequencing platform             : ${params.rgPlatform}"
        if ( params.rgLibrary ) log.info "Sequencing library              : ${params.rgLibrary}"
        if ( params.rgCenterName ) log.info "Sequencing center               : ${params.rgCenterName}"
        if ( params.rgDesc ) log.info "@RG Descritpiton                : ${params.rgDesc}"
        log.info ""
    }
    if ('bigwig' in params.stepList) {
        log.info "Bigwig parameters"
        log.info "-----------------"
        log.info "Tool                            : ${params.bigwigTool}"
        if ( params.bigwigTool.toLowerCase() == 'star' ) {
          def refPrefix = 'all'
          if ( 'wigRefPrefix' in paramNames ) {
            refPrefix = ['','-'].contains(params.wigRefPrefix) ? 'all' : params.wigRefPrefix
          }
          log.info "References prefix               : ${refPrefix}"
        }
        log.info ""
    }

    if ('quantification' in params.stepList) {
        log.info "Quantification parameters"
        log.info "-------------------------"
        log.info "Tool                            : ${params.quantificationTool}"
        if ( params.quantificationTool.toLowerCase() == 'rsem' ) {
          log.info "Skip CI                         : ${'rsemSkipCi' in paramNames ? params.rsemSkipCi : false}"
          log.info "Save RSEM model plot            : ${'rsemPlotModel' in paramNames ? params.rsemPlotModel : false}"
        }
        if ( 'quantificationMode' in paramNames ) {
          log.info "Mode                            : ${params.quantificationMode}"
        }
        log.info ""
    }

    def msg = "Output files db"
    log.info "=" * msg.size()
    log.info msg + " -> ${file(params.dbFile)}"
    log.info "=" * msg.size()
    log.info ""
}

/*
 * Given the input index file returns the number of unique samples,
 * the number of unique runs, and the lines of the index.
 * Params:
 * - tsvFile: a file object representing the TSV file
 * - printLog: a boolean to specify whether to print a summary of the TSV file (default: true)
 */
def readTsv(tsvFile, boolean printLog = true) {
    def (samples, ids, lines) = [[], [], []]
    tsvFile.eachLine { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        samples << sampleId
        ids << "${sampleId}${runId}"
        lines << line
    }
    nIds = ids.unique().size()
    nSamples = samples.unique().size()
    doMerge = nIds > nSamples
    if ( printLog ) {
        log.info "Dataset information"
        log.info "-------------------"
        log.info "Number of sequenced samples     : ${nSamples}"
        log.info "Number of sequencing runs       : ${nIds}"
        log.info "Merging                         : ${ doMerge ? 'by sample' : 'none' }"
        log.info ""
    }
    [ doMerge, lines ]
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

def parseJSON(data) {
    j = new JsonSlurper()
    j.parseText(data)
}