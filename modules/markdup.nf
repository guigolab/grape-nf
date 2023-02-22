process markdup {

    tag "${id.replace(':', '_')}-${params.markdupTool}-${params.markdupToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd)


    output:
    tuple val(id), val(sample), val(type), val(view), path("${prefix}.bam"), val(pairedEnd)

    script:
    cpus = task.cpus
    memory = (task.memory ?: 2.GB).toMega()
    prefix = "${bam.baseName}.markdup"

    command = "${task.process}/${params.markdupTool}${params.removeDuplicates ? '-remove' : ''}"
    template(command)
}