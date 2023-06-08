process inferExp {

    tag params.tag

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), stdout

    script:
    prefix = "${annotation.name.split('\\.', 2)[0]}"
    tagSuffix = "-${params.inferExpTool}-${params.inferExpToolVersion}"
    command = "${params.inferExpTool}"
    threshold = params.inferExpThreshold

    template(command)
}