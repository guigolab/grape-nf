process bamStats {

    tag "${id.replace(':', '_')}-${params.bamStatsTool}-${params.bamStatsToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd)
    tuple val(species), path(annotation)


    output:
    tuple val(id), val(sample), val(type), val(views), path('*.json'), val(pairedEnd)
 

    script:
    cpus = task.cpus
    type = "json"
    prefix = "${sample}"
    views = "BamStats"
    maxBuf = params.bamStatsMaxBuf
    logLevel = params.bamStatsLogLevel
    command = "${task.process}/${params.bamStatsTool}"

    template(command)
}