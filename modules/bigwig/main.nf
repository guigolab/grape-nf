process bigwig {

    tag params.tag

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)


    output:
    tuple val(sample), val(id), path('*.bw'), val(type), val(views), val(pairedEnd), val(readStrand)


    script:
    cpus = task.cpus
    type = "bigWig"
    prefix = "${sample}"
    tagSuffix = "-${params.bigwigTool}-${params.bigwigToolVersion}"
    wigRefPrefix = params.wigRefPrefix != "-" ? params.wigRefPrefix : ""
    views = params.bigwigViews[readStrand]
    command = "${params.bigwigTool}-${readStrand}"

    template(command)

}