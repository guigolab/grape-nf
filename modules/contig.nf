process contig {

    tag "${id.replace(':', '_')}-${params.contigTool}-${params.contigToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd), val(readStrand)
    tuple val(species), path(genomeFai)


    output:
    tuple val(id), val(sample), val(type), val(view), path('*.bed'), val(pairedEnd), val(readStrand)


    script:
    cpus = task.cpus
    type = 'bed'
    view = 'Contigs'
    prefix = "${sample}.contigs"
    command = "${task.process}/${params.contigTool}-${readStrand}"

    template(command)

}