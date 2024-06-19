params.rseqcVersion = '2.6.4--py27hf8a1672_2'
params.ucscVersion = '447--h2a80c09_1'
params.gtfToGenePredContainer = "${params.containerRepo}/ucsc-gtftogenepred:${params.ucscVersion}"
params.genePredContainer = "${params.containerRepo}/ucsc-genepredtobed:${params.ucscVersion}"
params.rseqcContainer = "public.ecr.aws/biocontainers/rseqc:${params.rseqcVersion}"
params.inferExpThreshold = '0.8'

process gtfToGenePred {

    tag "${prefix}"
    container params.gtfToGenePredContainer

    input:
    path(annotation)

    output:
    path("${prefix}.genePred")

    script:
    prefix = annotation.simpleName

    """
    gtfToGenePred ${annotation} -allErrors -ignoreGroupsWithoutExons ${prefix}.genePred 2> ${prefix}.genePred.err
    """
}

process genePredToBed {

    tag "${prefix}"
    container params.genePredContainer

    input:
    path(genePred)

    output:
    path("${prefix}.bed")

    script:
    prefix = "${genePred.simpleName}"

    """
    genePredToBed ${genePred} ${prefix}.bed
    """
}

process inferExp {

    tag "${sample}"
    container params.rseqcContainer

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), stdout

    script:
    prefix = "${annotation.simpleName}"

    """
    set -o pipefail

    grape_infer_experiment.py -i ${bam} -r ${annotation} --threshold ${params.inferExpThreshold} | tr -d '\\n'
    """
}
