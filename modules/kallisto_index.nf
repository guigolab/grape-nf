process kallisto_index {

    
    tag "${species}-kallisto"

    input:
    tuple val (species), path(kallisto_transcript_fa)



    output:
    tuple val(species), path{ "kallisto_index" }


    script:
    """
    #!/bin/bash
    kallisto index --index=kallisto_index whatever.transcripts.fa
    """

}