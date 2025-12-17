params.samtoolsVersion = '1.3.1--h0cf4675_11'
params.container = "quay.io/biocontainers/samtools:${params.samtoolsVersion}"

process sortBam {
    tag "${sample}"
    container params.container

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)

    script:
    memory = ((task.memory ?: 1.GB) / task.cpus) * 0.8
    prefix = "${bam.baseName}_sorted"

    """
    samtools sort -@ ${task.cpus} \
              -m ${memory.toMega()}M \
              -o ${prefix}.bam \
              ${bam}
    """
}
