process kallisto_transcript_fasta {

    label 'quantification'
    tag "${species}-kallisto"

    input:
    tuple val (species), path(kallisto_gtf)



    output:
    tuple val(species), path { "whatever.transcripts.fa" }


    script:
    """
    #!/bin/bash
    rsem-prepare-reference --gtf ${kallisto_gtf} /Users/pmarangio/PycharmProjects/grape-nf/data/genome.fa whatever
    """

}