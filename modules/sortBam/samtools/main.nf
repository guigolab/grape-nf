params.samtoolsVersion = '1.19.2--h50ea8bc_1'
params.container = "quay.io/biocontainers/samtools:${params.samtoolsVersion}"

process sortBam {
    tag "${sample}"
    container params.container

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)
    tuple val(sample), val(id), path("${prefix}.bam.bai"), val(type), val(view), val(pairedEnd)

    script:
    memory = (task.memory ?: 1.GB).toBytes() / task.cpus
    prefix = "${bam.baseName}_sorted"

    """
    samtools sort -@ ${task.cpus} \
              -m ${memory} \
              -o ${prefix}.bam \
              ${bam}
    samtools index ${prefix}.bam
    """
}

