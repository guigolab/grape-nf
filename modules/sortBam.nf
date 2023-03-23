process sortBam {
    tag "${id}-${params.mergeBamTool}-${params.mergeBamToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd)

    output:
    tuple val(id), val(sample), val(type), val(view), path("${prefix}.bam"), val(pairedEnd)

    script:
    cpus = task.cpus
    taskMemory = task.memory ?: 1.GB
    totalMemory = taskMemory.toBytes()
    threadMemory = totalMemory/cpus
    prefix = "${bam.baseName}_sorted"
    command = "${task.process}/${params.mergeBamTool}"

    template(command)
}

