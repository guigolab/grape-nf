def comprExts = ['gz', 'bz2', 'zip']

process index {

    label "mapping"
    tag "${params.mappingTool}-${params.mappingToolVersion}"

    input:
    path(genome)
    path(annotation)

    output:
    path("genomeDir")

    script:
    cpus = task.cpus
    memory = (task.memory ?: 1.GB).toBytes()
    sjOverHang = params.sjOverHang
    readLength = params.readLength
    genomeCompressed = genome.extension in comprExts ? "-genome-${genome.extension}" : ''
    annoCompressed = annotation.extension in comprExts ? "-anno-${annotation.extension}" : ''
    command = "${params.mappingTool}${genomeCompressed}${annoCompressed}"

    template(command)

}