params.samtoolsVersion = '1.19.2--h50ea8bc_1'
params.container = "${params.containerRepo}/samtools:${params.samtoolsVersion}"

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