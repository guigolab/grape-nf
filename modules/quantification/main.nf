process quantification {

    label 'quantification'
    tag params.tag

    input:
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)
    path(quantRef)


    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    cpus = task.cpus
    prefix = "${sample}"
    tagSuffix = "-${params.quantificationTool}-${params.quantificationToolVersion}"
    type = params.quantificationFileType
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    memory = (task.memory ?: 1.GB).toMega()
    command = "${params.quantificationTool}"
    if ( params.quantificationTool == 'RSEM') {
        command += "-${pairedEnd ? 'Paired-End' : 'Single-End'}"
    }
    command += "-${readStrand}"

    template(command)

}