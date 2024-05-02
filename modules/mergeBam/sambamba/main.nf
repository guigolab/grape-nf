params.sambambaVersion = '0.7.1--h984e79f_3'
params.container = "quay.io/biocontainers/sambamba:${params.sambambaVersion}"

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
    sambamba merge -t ${task.cpus} \\
               ${prefix}.bam \\
               *.bam && \\
    sambamba index ${prefix}.bam
    """

}
