process bigwig {

    tag "${id.replace(':', '_')}-${params.bigwigTool}-${params.bigwigToolVersion}"

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)


    output:
    tuple val(sample), val(id), path('*.bw'), val(type), val(views), val(pairedEnd), val(readStrand)


    script:
    cpus = task.cpus
    type = "bigWig"
    prefix = "${sample}"
    wigRefPrefix = params.wigRefPrefix != "-" ? params.wigRefPrefix : ""
    views = params.bigwigViews[readStrand]
    command = "${params.bigwigTool}-${readStrand}"

    template(command)

}