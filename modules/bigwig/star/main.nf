params.starVersion = '2.4.0j'
params.container = "grapenf/bigwig:star-${params.starVersion}"
params.wigRefPrefix = '-'

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
    
    def str1Prefix = null
    def str2Prefix = null
    def outWigStrand = null
    switch (readStrand) {
        case ~/(MATE1_)?SENSE/:
            str1Prefix = 'plusRaw'
            str2Prefix = 'minusRaw'
            outWigStrand = 'Stranded'
            views = ['MinusRawSignal', 'PlusRawSignal', 'MultipleMinusRawSignal','MultiplePlusRawSignal']
            break
        case ~/(MATE2_|ANTI)SENSE/:
            str1Prefix = 'minusRaw'
            str2Prefix = 'plusRaw'
            outWigStrand = 'Stranded'
            views = ['MinusRawSignal', 'PlusRawSignal', 'MultipleMinusRawSignal','MultiplePlusRawSignal']
            break
        default:
            str1Prefix = 'raw'
            str2Prefix = ''
            outWigStrand = 'Unstranded'
            views = ['RawSignal','MultipleRawSignal']
    }

    def cmd = []
    cmd << 'mkdir Signal'
    cmd << """\
        STAR --runThreadN ${task.cpus} \\
             --runMode inputAlignmentsFromBAM \\
             --inputBAMfile ${bam} \\
             --outWigType bedGraph \\
             --outWigStrand ${outWigStrand} \\
             --outWigReferencesPrefix ${params.wigRefPrefix} \\
             --outFileNamePrefix ./Signal/""".stripIndent()

    cmd << """\
        bedGraphToBigWig Signal/Signal.UniqueMultiple.str1.out.bg \\
                        ${genomeFai} \\
                        ${prefix}.UniqueMultiple.${str1Prefix}.bw
        bedGraphToBigWig Signal/Signal.Unique.str1.out.bg \\
                        ${genomeFai} \\
                        ${prefix}.Unique.${str1Prefix}.bw""".stripIndent()
    if ( outWigStrand == 'Stranded' ) {
        cmd << """\
            bedGraphToBigWig Signal/Signal.UniqueMultiple.str2.out.bg \\
                            ${genomeFai} \\
                            ${prefix}.UniqueMultiple.${str2Prefix}.bw
            bedGraphToBigWig Signal/Signal.Unique.str2.out.bg \\
                            ${genomeFai} \\
                            ${prefix}.Unique.${str2Prefix}.bw""".stripIndent()
    }
    cmd.join('\n')
}