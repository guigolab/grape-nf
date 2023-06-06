process inferExp {

    tag "${id.replace(':', '_')}-${params.inferExpTool}-${params.inferExpToolVersion}"

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), stdout

    script:
    prefix = "${annotation.name.split('\\.', 2)[0]}"
    command = "${params.inferExpTool}"
    threshold = params.inferExpThreshold

    template(command)
}