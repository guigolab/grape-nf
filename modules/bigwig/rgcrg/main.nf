params.rgcrgVersion = '0.1'
params.container = "grapenf/bigwig:rgcrg-${params.rgcrgVersion}"

process bigwig {

    tag "${sample}"
    container params.container

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path('*.bw'), val(type), val(views), val(pairedEnd), val(readStrand)

    script:
    type = "bigWig"
    prefix = "${sample}"
    views = ['RawSignal']

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
        cmd << "genomeCoverageBed -split -bg -ibam ${bam} > ${prefix}.${str1Prefix}.bedgraph"
    } else {
        cmd << """\
            sambamba view -h -t ${task.cpus} ${bam} \\
            | awk 'BEGIN {OFS="\\t"}${awkStr}' \\
            | sambamba view -S -f bam -l 0 /dev/stdin \\
            | tee >(
                genomeCoverageBed -strand + -split -bg -ibam - \\
                > ${prefix}.${str1Prefix}.bedgraph
                ) \\
            | genomeCoverageBed -strand - -split -bg -ibam - \\
            > ${prefix}.${str2Prefix}.bedgraph""".stripIndent()
    }

    cmd << """\
        bedGraphToBigWig ${prefix}.${str1Prefix}.bedgraph ${genomeFai} ${prefix}.${str1Prefix}.bw \
    """.stripIndent()
    if ( str2Prefix ) {
        cmd << """\
            bedGraphToBigWig ${prefix}.${str2Prefix}.bedgraph ${genomeFai} ${prefix}.${str2Prefix}.bw \
        """.stripIndent()
    }

    cmd.join('\n')
}