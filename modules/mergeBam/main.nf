def pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

process mergeBam {

    tag "${id.replace(':', '_')}-${params.mergeBamTool}-${params.mergeBamToolVersion}"

    input:
    tuple val(sample), val(id), path("${sample}_??.bam"), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), path("${prefix}.bam"), val(type), val(view), val(pairedEnd)

    script:
    cpus = task.cpus
    id = id.sort().join(':')
    prefix = "${sample}${pref}_to${view.replace('Alignments','')}"
    command = "${params.mergeBamTool}"

    template(command)

}