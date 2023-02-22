def comprExts = ['gz', 'bz2', 'zip']

process txIndex {

    label 'quantification'
    tag "${species}-${params.quantificationTool}-${params.quantificationToolVersion}"

    input:
    tuple val(species), path(genome)
    tuple val(species), path(annotation)

    output:
    tuple val(species), path('txDir')

    script:
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${task.process}/${params.quantificationTool}${genomeCompressed}${annoCompressed}"

    template(command)

}