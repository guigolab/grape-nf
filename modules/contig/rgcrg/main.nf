params.rgcrgVersion = '0.1'
params.container = "grapenf/contig:rgcrg-${params.rgcrgVersion}"

process contig {

    tag "${sample}"
    container params.container

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)


    output:
    tuple val(sample), val(id), path('*.bed'), val(type), val(view), val(pairedEnd), val(readStrand)


    script:
    type = 'bed'
    view = 'Contigs'
    prefix = "${sample}.contigs"

    def str1Prefix = 'raw'
    def str2Prefix = ''
    def awkStr = ''

    if ( readStrand in  [ 'SENSE', 'ANTISENSE', 'MATE1_SENSE', 'MATE2_SENSE' ] ) {
        str1Prefix = 'plusRaw'
        str2Prefix = 'minusRaw'
        views = ['MinusRawSignal', 'PlusRawSignal']
    }
    switch (readStrand) {
        case 'ANTISENSE':
            awkStr = '{if ($1!~/^@/) {$2=xor($2,0x10)}; print}'
            break
        case 'MATE1_SENSE':
            awkStr = '{if ($1!~/^@/ && and($2,128)>0) {$2=xor($2,0x10)}; print}'
            break
        case 'MATE2_SENSE':
            awkStr = '{if ($1!~/^@/ && and($2,64)>0) {$2=xor($2,0x10)}; print}'
            break            
    }

    def cmd = []

    if ( ! str2Prefix ) {
        cmd << """\
            sambamba view  -f bam -F '[NH] == 1' ${bam} \\
            | bamToBed -i - \\
            | sort -T. -k1,1 -k2,2n \\
            | mergeBed \\
            > ${prefix}.bed""".stripIndent()
    } else {
        cmd << """\
            sambamba view -h -t ${task.cpus} ${bam} \\
            | awk 'BEGIN {OFS="\\t"}${awkStr}' \\
            | sambamba view -S -f bam -l 0 -F '[NH] == 1' /dev/stdin \\
            | tee >(
                genomeCoverageBed -strand + -split -bg -ibam - \\
                > ${prefix}.${str1Prefix}.bedgraph
                ) \\
            | genomeCoverageBed -strand - -split -bg -ibam - \\
            > ${prefix}.${str2Prefix}.bedgraph""".stripIndent()
        cmd << """\
            contigsNew.py --chrFile ${genomeFai} \
                        --fileP ${prefix}.plusRaw.bedgraph \
                        --fileM ${prefix}.minusRaw.bedgraph \
                        --sortOut \
            > ${prefix}.bed""".stripIndent()
    }


    cmd.join('\n')
}