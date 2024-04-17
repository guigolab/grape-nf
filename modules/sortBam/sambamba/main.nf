params.sambambaVersion = '0.7.1'
params.container = "grapenf/sambamba:${params.sambambaVersion}"

process sortBam {
    tag "${sample}"
    container params.container

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)
    tuple val(sample), val(id), path("${prefix}.bam.bai"), val(type), val(view), val(pairedEnd)

    script:
    memory = (task.memory ?: 1.GB).toBytes()
    prefix = "${bam.baseName}_sorted"

    """
    sambamba sort -t ${task.cpus} \
              -m ${memory} \
              -o ${prefix}.bam \
              ${bam}
    sambamba index ${prefix}.bam
    """
}

