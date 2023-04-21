process kallisto_chrom_length {

    tag "${species}-kallisto"

    input:
    tuple val(species), path(genome)


    output:
    tuple val(species), path{ "${genome.name.replace('.fa.fai','')}.size" }


    script:
    """
    #!/bin/bash
    cut -f1,2 ${genome} > ${genome.name.replace('.fa.fai','')}.size
    """

}