params.rsemVersion = '1.2.21'
params.container = "grapenf/quantification:rsem-${params.rsemVersion}"
params.rsemSkipCi = false
params.rsemPlotModel = false

process index {

    tag "${genome}-${annotation}"
    container params.container

    input:
    path(genome)
    path(annotation)

    output:
    path('txDir')

    script:
    def cmd = []

    cmd << "mkdir txDir"
    if ( genome.extension in params.comprExts ) {
        cmd << """\
            mkfifo ${genome.baseName}
            zcat ${genome} > ${genome.baseName} &
        """.stripIndent()
        genome = ${genome.baseName}
    }
    if ( annotation.extension in params.comprExts ) {
        cmd << """\
            mkfifo ${annotation.baseName}
            zcat ${annotation} > ${annotation.baseName} &
        """.stripIndent()
        annotation = ${annotation.baseName}
    }
    cmd << """\
        rsem-prepare-reference --gtf ${annotation} \\
                               ${genome} \\
                               txDir/RSEMref""".stripIndent()
    cmd.join('\n')
}

process quantify {

    tag "${sample}"
    container params.container

    input:
    path(quantRef)
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    prefix = "${sample}"
    type = 'tsv'
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    def memory = (task.memory ?: 1.GB).toMega()
    def forwardProb = null

    switch (readStrand) {
        case ~/(ANTI|MATE2_)SENSE/:
            forwardProb = '0'
            break
        case ~/(MATE1_)?SENSE/:
            forwardProb = '1'
            break
    }
    
    def cmd = []
    cmd << "sambamba sort -t ${task.cpus} -m ${memory}MB -N -M -l 0 -o - ${bam} \\"
    cmd << """\
        | rsem-calculate-expression -p ${task.cpus} \\
                                    --bam \\
                                    --seed 12345 \\
                                    --estimate-rspd  \\
                                    --no-bam-output \\""".stripIndent()
    if ( pairedEnd ) {
        cmd << """\
                            --paired-end \\"""
    }
    if ( forwardProb ) {
        cmd << """\
                            --forward-prob ${forwardProb} \\"""
    }
    if ( ! params.rsemSkipCi ) {
        cmd << """\
                            --calc-ci \\
                            --ci-memory ${memory} \\"""
    }
    cmd << """\
                            - \\
                            ${quantRef}/RSEMref \\
                            ${prefix}"""
    if ( params.rsemPlotModel ) {
        cmd << "rsem-plot-model ${prefix} ${prefix}.pdf"
    }
    cmd.join('\n')
}