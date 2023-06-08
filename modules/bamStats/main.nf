process bamStats {

    tag params.tag

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)


    output:
    tuple val(sample), val(id), path('*.json'), val(type), val(views), val(pairedEnd)
 

    script:
    cpus = task.cpus
    type = "json"
    prefix = "${sample}"
    tagSuffix = "-${params.bamStatsTool}-${params.bamStatsToolVersion}"
    views = "BamStats"
    maxBuf = params.bamStatsMaxBuf
    logLevel = params.bamStatsLogLevel
    command = "${params.bamStatsTool}"

    template(command)
}