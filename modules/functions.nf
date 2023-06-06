// imports
import groovy.json.JsonSlurper

// Print pipeline usage and exit if `--help` is passed
def printUsage() {
    if (params.help) {
        log.info ''
        log.info 'G R A P E ~ RNA-seq Pipeline'
        log.info '----------------------------'
        log.info 'Run the GRAPE RNA-seq pipeline on a set of data.'
        log.info ''
        log.info 'Usage: '
        log.info '    grape-pipeline --index INDEX_FILE --genome GENOME_FILE --annotation ANNOTATION_FILE [OPTION]...'
        log.info ''
        log.info 'Options:'
        log.info '    --help                              Show this message and exit.'
        log.info '    --index INDEX_FILE                  Index file.'
        log.info '    --genome GENOME_FILE                Reference genome file(s).'
        log.info '    --annotation ANNOTAION_FILE         Reference gene annotation file(s).'
        log.info '    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all'
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
        log.info ''
        exit 1
    }
}

// Print pipeline log
def printLog() {
    log.info ""
    log.info "G R A P E ~ RNA-seq Pipeline"
    log.info ""
    log.info "General parameters"
    log.info "------------------"
    log.info "Index file                      : ${params.index}"
    log.info "Genome                          : ${params.genome}"
    log.info "Annotation                      : ${params.annotation}"
    log.info "Pipeline steps                  : ${params.stepList.join(" ")}"
    log.info "Pipeline profile                : ${workflow.profile}"
    log.info ""

    if ('mapping' in params.stepList) {
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
    if ('bigwig' in params.stepList) {
        log.info "Bigwig parameters"
        log.info "-----------------"
        // log.info "Tool                            : ${bigwigTool}"
        log.info "References prefix               : ${['','-'].contains(params.wigRefPrefix) ? 'all' : params.wigRefPrefix}"
        log.info ""
    }

    if ('quantification' in params.stepList) {
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
    log.info msg + " -> ${file(params.dbFile)}"
    log.info "=" * msg.size()
    log.info ""
}

/*
 * Given the input index file returns the number of unique samples,
 * the number of unique runs, and the lines of the index.
 * Params:
 * - tsvFile: a file object representing the TSV file
 */
def readTsv(tsvFile, boolean printLog = true) {
    def (samples, ids, lines) = [[], [], []]
    tsvFile.eachLine { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        samples << sampleId
        ids << runId // TODO: fix wrong number of runs
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