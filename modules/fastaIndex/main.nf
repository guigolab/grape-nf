def comprExts = ['gz', 'bz2', 'zip']

process fastaIndex {

    tag "${params.fastaIndexTool}-${params.fastaIndexToolVersion}"

    input:
    path(genome)

    output:
    path { "${genome.name.replace('.gz','')}.fai" }

    script:
    compressed = genome.extension in comprExts ? "-${genome.extension}" : ''
    command = "${params.fastaIndexTool}${compressed}"
    
    template(command)

}