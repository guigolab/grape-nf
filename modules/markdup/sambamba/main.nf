params.sambambaVersion = '0.7.1--h984e79f_3'
params.container = "quay.io/biocontainers/sambamba:${params.sambambaVersion}"

process markdup {

    tag "${sample}"
    container params.container

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)


    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd), emit: dedupedAlignments
    tuple val(sample), val(id), path("${prefix}.bam.bai"), val('bai'), val("${view}Index"), val(pairedEnd), emit: dedupedAlignmentsIndices

    script:
    memory = (task.memory ?: 2.GB).toMega()
    prefix = "${bam.baseName}.markdup"

    def cmd = []
    cmd << """\
        sambamba markdup -t ${task.cpus} \\
                         --sort-buffer-size ${memory} \\""".stripIndent()
    if ( params.removeDuplicates ) {
        cmd << '                 --remove-duplicates \\'
    }
    cmd << """\
                         ${bam} \\
                         ${prefix}.bam && \\
        sambamba index ${prefix}.bam""".stripIndent()

    cmd.join('\n')
}
