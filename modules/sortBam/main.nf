process sortBam {
    tag params.tag

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)

    script:
    cpus = task.cpus
    taskMemory = task.memory ?: 1.GB
    totalMemory = taskMemory.toBytes()
    threadMemory = totalMemory/cpus
    prefix = "${bam.baseName}_sorted"
    tagSuffix = "-${params.mergeBamTool}-${params.mergeBamToolVersion}"
    command = "${params.mergeBamTool}"

    template(command)
}

