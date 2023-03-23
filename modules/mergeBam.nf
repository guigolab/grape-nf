def pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

process mergeBam {

    tag "${id.replace(':', '_')}-${params.mergeBamTool}-${params.mergeBamToolVersion}"

    input:
    tuple val(id), val(sample), val(type), val(view), path("${sample}_??.bam"), val(pairedEnd)

    output:
    tuple val(id), val(sample), val(type), val(view), path("${prefix}.bam"), val(pairedEnd)

    script:
    cpus = task.cpus
    id = id.sort().join(':')
    prefix = "${sample}${pref}_to${view.replace('Alignments','')}"
    command = "${task.process}/${params.mergeBamTool}"

    template(command)

}