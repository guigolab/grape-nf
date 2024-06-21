include { parseJSON } from "../../functions"

params.samtoolsVersion = '1.19.2--h50ea8bc_1'
params.samtoolsContainer = "${params.containerRepo}/samtools:${params.samtoolsVersion}"
params.bamstatsVersion = '0.3.5--he881be0_0'
params.bamStatsContainer = "quay.io/biocontainers/bamstats:${params.bamstatsVersion}"
params.bamStatsMaxBuf = '1000000'
params.bamStatsLogLevel = 'info'

process getProtocol {

    tag "${sample}"
    container params.bamStatsContainer

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view)


    output:
    tuple val(sample), val(id),  stdout
 

    script:
    prefix = "${sample}"

    """
    bamstats -c ${task.cpus} \\
             -i ${bam} \\
             -n 1000 \\
             --max-buf ${params.bamStatsMaxBuf} \\
             --loglevel ${params.bamStatsLogLevel}
    """
}

process toFastq {

    tag "${sample}"
    container params.samtoolsContainer

    input:
    tuple val(sample), val(id), path(bam), val(type), val(view), val(protocol)

    output:
    tuple val(sample), val(id), path("${prefix}*.fastq.gz"), val(type), val(view)

    script:
    prefix = "${id}"
    type = "fastq"
    view = "FqRd"

    def cmd = []
    def inputBam = bam
    outParams = "-0 ${prefix}.fastq.gz"
    if ( protocol == "PairedEnd" ) {
	view = ['FqRd1', 'FqRd2']
	inputBam = "${prefix}.collate.bam"
	outParams = "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz"
	cmd << "mkfifo ${prefix}.collate.bam"
	cmd << "samtools collate --threads ${task.cpus} ${bam} ${prefix}.collate &"
    }
     
    cmd << "samtools fastq -N -c 9 --threads ${task.cpus} ${outParams} ${inputBam}"

    cmd.join('\n')
}

workflow bamToFastq {
    take:
      genomeAlignments

    main:
    
    getProtocol( genomeAlignments )
    getProtocol.out.map {
	d = parseJSON(it[-1])
	it[0..1] + [ d.general.protocol ]
    }.set { protocol }

    genomeAlignments
	.join( protocol, by: [0,1] )
	.set { bams }
    
    toFastq( bams )

    emit:
      toFastq.out
}
