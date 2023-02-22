process inferExp {

    tag "${id.replace(':', '_')}-${params.inferExpTool}-${params.inferExpToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd)
    tuple val(species), path(annotation)

    output:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd), stdout
    

    script:
    prefix = "${annotation.name.split('\\.', 2)[0]}"
    command = "${task.process}/${params.inferExpTool}"
    threshold = params.inferExpThreshold

    template(command)
}