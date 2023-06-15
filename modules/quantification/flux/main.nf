params.fluxVersion = '1.6.1'
params.container = "grapenf/quantification:flux-${params.fluxVersion}"
params.fluxProfile = false

process quantify {

    tag "${sample}"
    container params.container

    input:
    path(annotation)
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    prefix = "${sample}"
    type = 'gtf'
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    def memory = (task.memory ?: 1.GB).toMega()
    def mode = pairedEnd ? 'PAIRED' : 'SINGLE'

    def cmd = []
    cmd << "samtools index ${bam}"
    if ( params.fluxProfile ) {
        cmd << """\
            flux-capacitor --profile \\
                        -i ${bam} \\ 
                        -a ${annotation} \\""".stripIndent()
        if ( readStrand in ['SENSE', 'ANTISENSE', 'MATE1_SENSE', 'MATE2_SENSE'] ) {
            cmd << """\
            -m ${mode}_STRANDED \\
            --read-strand ${readStrand} \\"""
        }
        cmd << "            --profile-file ${prefix}_profile.json"
    }
    cmd << """\
        flux-capacitor -i ${bam} \\
                    -a ${annotation} \\""".stripIndent()
    if ( readStrand in ['SENSE', 'ANTISENSE', 'MATE1_SENSE', 'MATE2_SENSE'] ) {
        cmd << """\
            -m ${mode}_STRANDED \\
            --read-strand ${readStrand} \\"""
    }
    cmd << "            -o ${prefix}.isoforms.gtf"
    cmd << """\
        TrtoGn_RPKM.sh -a ${annotation} \\
                    -i ${prefix}.isoforms.gtf \\
                    -o ${prefix}.genes.gff""".stripIndent()
    cmd.join('\n')
}