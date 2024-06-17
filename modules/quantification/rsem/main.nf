params.rsemVersion = '1.3.3--pl5321h6b7c446_6'
params.container = "${params.containerRepo}/rsem:${params.rsemVersion}"
params.rsemCalcCI = false
params.rsemPlotModel = false

process index {

    tag "${genome.simpleName}-${annotation.simpleName}"
    container params.container

    input:
    path(genome)
    path(annotation)

    output:
    path('txDir')

    script:
    def cmd = []
    def genomeFile = genome.name
    def annotationFile = annotation.name

    cmd << "mkdir txDir"
    if ( genome.extension in params.comprExts ) {
        genomeFile = genome.baseName
        cmd << """\
            mkfifo ${genomeFile}
            zcat ${genome} > ${genomeFile} &
        """.stripIndent()
    }
    if ( annotation.extension in params.comprExts ) {
        annotationFile = annotation.baseName
        cmd << """\
            mkfifo ${annotationFile}
            zcat ${annotation} > ${annotationFile} &
        """.stripIndent()
    }
    cmd << """\
        rsem-prepare-reference --gtf ${annotationFile} \\
                               ${genomeFile} \\
                               txDir/RSEMref""".stripIndent()
    cmd.join('\n')
}

process validateBam {
    tag "${sample}"
    container params.container

    input:
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), stdout

    """
    rsem-sam-validator ${bam}
    """
}

process calculateExpression {

    tag "${sample}"
    container params.container

    input:
    path(quantRef)
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand), val(convertBam)

    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    prefix = "${sample}"
    rsemPrefix = bam.simpleName
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
    if ( convertBam ) {
        rsemPrefix = "${rsemPrefix}.rsem"
        def convertOpts = ""
        if ( params.rsemVersion.startsWith("1.3") ) {
            def convertMemory = Math.max(1024, ( memory / task.cpus / 2 ) as long)
            convertOpts = "-p ${task.cpus} --memory-per-thread ${convertMemory}M"
        }
        cmd << """\
            mkfifo ${rsemPrefix}.bam
            convert-sam-for-rsem ${convertOpts} \\
                                 ${bam} \\
                                 ${rsemPrefix} &
        """.stripIndent()
    }
    cmd << """\
        rsem-calculate-expression -p ${task.cpus} \\
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
    if ( params.rsemCalcCI ) {
        def ciMemory = Math.max(1024, (memory / 4) as long)
        cmd << """\
                            --calc-ci \\
                            --ci-memory ${ciMemory} \\"""
    }
    cmd << """\
                            ${rsemPrefix}.bam \\
                            ${quantRef}/RSEMref \\
                            ${prefix}"""
    if ( params.rsemPlotModel ) {
        cmd << "rsem-plot-model ${prefix} ${prefix}.pdf"
    }
    cmd.join('\n')
}
