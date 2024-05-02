params.samtoolsVersion = '1.19.2--h50ea8bc_1'
params.container = "quay.io/biocontainers/samtools:${params.samtoolsVersion}"

process mergeBam {

    tag "${sample}-${id}"
    container params.container

    input:
    tuple val(sample), val(id), path("${sample}_??.bam"), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)

    script:
    id = id.sort().join(':')
    prefix = "${sample}_m${params.maxMismatches}_n${params.maxMultimaps}_to${view.replace('Alignments','')}"

    """
    samtools merge -@ ${task.cpus} \\
                -p \\
                ${prefix}.bam \\
                *.bam && \\
    samtools index ${prefix}.bam
    """

}
