def comprExts = ['gz', 'bz2', 'zip']

process txIndex {

    label 'quantification'
    tag "${params.quantificationTool}-${params.quantificationToolVersion}"

    input:
    path(genome)
    path(annotation)

    output:
    path('txDir')

    script:
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${params.quantificationTool}${genomeCompressed}${annoCompressed}"

    template(command)

}