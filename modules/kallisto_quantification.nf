process kallisto_quantification {


    tag "${species}-kallisto"

    input:

    tuple val (species), path(index)
    tuple val (species), path(chromosomes)
    tuple val (species), path(gtf)
    tuple path(readone), path(readtwo)



    output:
    tuple val(species), path("outputdir")



    script:
    """
    #!/bin/bash
    
    kallisto quant -i ${index} -c ${chromosomes} --gtf ${gtf} -o outputdir ${readone} ${readtwo}
    """

}