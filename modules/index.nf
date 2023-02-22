def comprExts = ['gz', 'bz2', 'zip']

process index {

    label "mapping"
    tag "${species}-${params.mappingTool}-${params.mappingToolVersion}"

    input:
    tuple val(species), path(genome)
    tuple val(species), path(annotation)

    output:
    tuple val(species), path("genomeDir")

    script:
    cpus = task.cpus
    memory = (task.memory ?: 1.GB).toBytes()
    sjOverHang = params.sjOverHang
    readLength = params.readLength
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${task.process}/${params.mappingTool}${genomeCompressed}${annoCompressed}"

    template(command)

}