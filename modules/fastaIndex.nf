def comprExts = ['gz', 'bz2', 'zip']

process fastaIndex {

    tag "${species}-${params.fastaIndexTool}-${params.fastaIndexToolVersion}"

    input:
    tuple val (species), path(genome)
    tuple val (species), path(annotation)

    output:
    tuple val(species), path { "${genome.name.replace('.gz','')}.fai" }

    script:
    compressed = genome.extension in comprExts ? "-${genome.extension}" : ''
    command = "${task.process}/${params.fastaIndexTool}${compressed}"
    template(command)

}