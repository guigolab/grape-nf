process quantification {

    label 'quantification'
    tag "${id.replace(':', '_')}-${params.quantificationTool}-${params.quantificationToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd), val(readStrand)
    tuple val(species), val(quantRef)


    output:
    tuple val(id), val(sample), val(type), val(viewTx), path("*isoforms*"), val(pairedEnd), val(readStrand)
    tuple val(id), val(sample), val(type), val(viewGn), path("*genes*"), val(pairedEnd), val(readStrand)



    script:
    cpus = task.cpus
    prefix = "${sample}"
    type = params.quantificationFileType
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    memory = (task.memory ?: 1.GB).toMega()
    command = "${task.process}/${params.quantificationTool}"
    if ( params.quantificationTool == 'RSEM') {
        command += "-${pairedEnd ? 'Paired-End' : 'Single-End'}"
    }
    command += "-${readStrand}"

    template(command)

}