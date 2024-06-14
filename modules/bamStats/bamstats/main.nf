params.bamstatsVersion = '0.3.5--he881be0_0'
params.container = "${params.containerRepo}/bamstats:${params.bamstatsVersion}"
params.bamStatsMaxBuf = '1000000'
params.bamStatsLogLevel = 'info'

process bamStats {

    tag "${sample}"
    container params.container

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)


    output:
    tuple val(sample), val(id), path('*.json'), val(type), val(view), val(pairedEnd)
 

    script:
    type = "json"
    view = "BamStats"
    prefix = "${sample}"

    """
    bamstats -c ${task.cpus} \\
             -i ${bam} \\
             -a ${annotation} \\
             -o ${prefix}_stats.json \\
             --max-buf ${params.bamStatsMaxBuf} \\
             --loglevel ${params.bamStatsLogLevel} \\
             -u
    """
}