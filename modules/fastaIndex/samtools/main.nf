params.samtoolsVersion = '1.3.1--h0cf4675_11'
params.container = "quay.io/biocontainers/samtools:${params.samtoolsVersion}"

process fastaIndex {

    tag "${genome.simpleName}"
    container params.container

    input:
    path(genome)

    output:
    path('*.fai')

    script:
    def compressed = genome.extension in params.comprExts
    def genomeFile = genome.name
    
    def cmd = []
    if ( compressed ) { 
        genomeFile = genome.baseName
        cmd << "gzip -dc ${genome} > ${genomeFile}"
    }
    cmd << "samtools faidx ${genomeFile}"
    if ( compressed ) {
        cmd << "rm ${genomeFile}"
    }
    cmd.join('\n')
}