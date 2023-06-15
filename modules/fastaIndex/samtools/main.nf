params.samtoolsVersion = '1.3.1--h0cf4675_11'
params.container = "quay.io/biocontainers/samtools:${params.samtoolsVersion}"

process fastaIndex {

    tag "${genome}"
    container params.container

    input:
    path(genome)

    output:
    path('*.fai')

    script:
    def compressed = genome.extension in params.comprExts
    
    def cmd = []
    if ( compressed ) { 
        cmd << "gzip -dc ${genome} > ${genome.baseName}"
        genome = genome.baseName
    }
    cmd << "samtools faidx ${genome}"
    if ( compressed ) {
        cmd << "rm ${genome}"
    }
    cmd.join('\n')
}