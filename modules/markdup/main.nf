process markdup {

    tag params.tag

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)


    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd), emit: dedupedAlignments
    tuple val(sample), val(id), path("${prefix}.bam.bai"), val('bai'), val("${view}Index"), val(pairedEnd), emit: dedupedAlignmentsIndices

    script:
    cpus = task.cpus
    memory = (task.memory ?: 2.GB).toMega()
    prefix = "${bam.baseName}.markdup"
    tagSuffix = "-${params.markdupTool}-${params.markdupToolVersion}"
    command = "${params.markdupTool}${params.removeDuplicates ? '-remove' : ''}"
    
    template(command)
}