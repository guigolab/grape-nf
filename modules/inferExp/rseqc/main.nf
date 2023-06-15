params.rseqcVersion = '2.6.4'
params.container = "grapenf/inferexp:rseqc-${params.rseqcVersion}"
params.inferExpThreshold = '0.8'

process inferExp {

    tag "${sample}"
    container params.container

    input:
    path(annotation)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd)

    output:
    tuple val(sample), val(id), stdout

    script:
    prefix = "${annotation.simpleName}"

    """
    set -o pipefail

    gtfToGenePred ${annotation} -allErrors -ignoreGroupsWithoutExons ${prefix}.genePred 2> ${prefix}.genePred.err
    genePredToBed ${prefix}.genePred ${prefix}.bed
    grape_infer_experiment.py -i ${bam} -r ${prefix}.bed --threshold ${params.inferExpThreshold} | tr -d '\\n'
    """
}