process contig {

    tag "${id.replace(':', '_')}-${params.contigTool}-${params.contigToolVersion}"

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)


    output:
    tuple val(sample), val(id), path('*.bed'), val(type), val(view), val(pairedEnd), val(readStrand)


    script:
    cpus = task.cpus
    type = 'bed'
    view = 'Contigs'
    prefix = "${sample}.contigs"
    command = "${params.contigTool}-${readStrand}"

    template(command)

}