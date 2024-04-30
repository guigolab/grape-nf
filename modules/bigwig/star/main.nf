params.starVersion = "2.4.0j--h9ee0642_2"
// params.starVersion = "2.7.10a--h9ee0642_0"
params.ucscVersion = '455--h2a80c09_1'
params.starContainer = "quay.io/biocontainers/star:${params.starVersion}"
params.bgtobwContainer = "quay.io/biocontainers/ucsc-bedgraphtobigwig:${params.ucscVersion}"
params.wigRefPrefix = '-'

process signal {
    tag "${sample}"
    container params.starContainer

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path('Signal'), val(type), val('SignalBedgraphs'), val(pairedEnd), val(readStrand)

    script:
    type = "bedGraph"
    outWigStrand = readStrand == "NONE" ? 'Unstranded' : 'Stranded'
    """
    mkdir Signal
    STAR --runThreadN ${task.cpus} \\
         --runMode inputAlignmentsFromBAM \\
         --inputBAMfile ${bam} \\
         --outWigType bedGraph \\
         --outWigStrand ${outWigStrand} \\
         --outWigReferencesPrefix ${params.wigRefPrefix} \\
         --outFileNamePrefix ./Signal/
    """
}

process bw {

    tag "${sample}"
    container params.bgtobwContainer

    input:
    path(genomeFai)
    tuple val(sample), val(id), path(signal), val(type), val(view), val(pairedEnd), val(readStrand)

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

workflow bigwig {
    take:
      fastaIndex
      genomeAlignments

    main:
      signal( fastaIndex, genomeAlignments )
      bw( fastaIndex, signal.out )

    emit:
      bw.out
}
