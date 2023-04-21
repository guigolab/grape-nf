process kallisto_prep_gtf {

    tag "${species}-kallisto"

    input:
    tuple val (species), path(annotation)
    //tuple val (species), path(genome)
    //set species, file(genome) from fastaIndexGenomes
    //set species, file(annotation) from fastaIndexAnnotations

    output:
    tuple val(species), path { "kallisto.gtf" }
    //set species, file { "${genome.name.replace('.gz','')}.fai" } into fastaIndexOutput

    script:
    """
    #!/bin/bash
    gawk '( \$3 ~ /gene/ )' ${annotation} > kallisto.gtf
    gawk '( \$3 ~ /transcript/ )' ${annotation}  >> kallisto.gtf
    gawk '( \$3 ~ /exon/ && \$7 ~ /+/ )' ${annotation} | sort -k1,1 -k4,4n >> kallisto.gtf
    gawk '( \$3 ~ /exon/ && \$7 ~ /-/ )' ${annotation}  | sort -k1,1 -k4,4nr >> kallisto.gtf
    """
}