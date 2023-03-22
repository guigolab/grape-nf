process bigwig {

    tag "${id.replace(':', '_')}-${params.bigwigTool}-${params.bigwigToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path(bam), val(pairedEnd), val(readStrand)
    tuple val(species), path(genomeFai)


    output:
    tuple val(id), val(sample), val(type), val(views), path('*.bw'), val(pairedEnd), val(readStrand)


    script:
    cpus = task.cpus
    type = "bigWig"
    prefix = "${sample}"
    wigRefPrefix = params.wigRefPrefix != "-" ? params.wigRefPrefix : ""
    views = params.bigwigViews[readStrand]
    command = "${task.process}/${params.bigwigTool}-${readStrand}"

    template(command)

}