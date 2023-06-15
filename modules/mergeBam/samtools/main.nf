params.samtoolsVersion = '1.3.1--h0cf4675_11'
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
    (
        samtools view -H *.bam | grep -v "@RG";
        for f in *.bam; do 
            samtools view -H \$f | grep "@RG";
        done
    ) > header.txt && \\
    samtools merge -@ ${task.cpus} \\
                -h header.txt \\
                ${prefix}.bam \\
                *.bam && \\
    samtools index ${prefix}.bam
    """

}