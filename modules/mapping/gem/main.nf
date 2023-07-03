params.gemtoolsVersion = "1.7.1"
params.container = "grapenf/mapping:gem-${params.gemtoolsVersion}"
params.mappingSortTool = 'samtools'
params.readLength = 150

process index {

    container params.container
    tag "${genome.simpleName}-${annotation.simpleName}"

    input:
    path(genome)
    path(annotation)

    output:
    path("genomeDir")

    script:
    def pigzCpus = Math.min(task.cpus, 4)
    def memory = (task.memory ?: 1.GB).toBytes()
    def compressedGenome = genome.extension in params.comprExts
    def compressedAnnotation = annotation.extension in params.comprExts
    def genomeFile = genome.name
    def annotationFile = annotation.name

    def cmd = [ 'mkdir genomeDir' ]
    if ( compressedGenome ) {
        genomeFile = genome.baseName
        cmd << "pigz -p ${pigzCpus} -dc ${genome} | sed 's/ .*//' > ${genomeFile}"
    }
    if ( compressedAnnotation ) {
        annotationFile = annotation.baseName
        cmd << "mkfifo ${annotationFile}"
        cmd << "pigz -p ${pigzCpus} -dc ${annotation} > ${annotationFile} &"
    }
    cmd << """\
        gemtools index -i ${genomeFile} \\
                       -t ${task.cpus} \\
                       -o genomeDir/genome_index.gem""".stripIndent()
    cmd << """\
        gemtools t-index -i genomeDir/genome_index.gem \\
                         -a ${annotationFile} \\
                         -t ${task.cpus} \\
                         -o genomeDir/transcript_index \\
                         -m ${params.readLength}""".stripIndent()
    if ( compressedGenome ) {
        cmd << "rm ${genomeFile}"
    }
    if ( compressedAnnotation ) {
        cmd << "rm ${annotationFile}"
    }
    cmd.join('\n')

}

process map {

    container params.container
    tag "${sample}-${id}"

    input:
    path(annotation)
    path(genomeDir)
    tuple val(sample), val(id), path(reads), val(type), val(view), val(qualityOffset)

    output:
    tuple val(sample), val(id), path("*toGenome.bam"), val(type), val("Genome${view}"), val(pairedEnd), emit: genomeAlignments
    tuple val(sample), val(id), path("*toTranscriptome.bam"), val(type), val("Transcriptome${view}"), val(pairedEnd), optional: true, emit: transcriptomeAlignments
    tuple val(sample), val(id), path("*toGenome.bam.bai"), val('bai'), val("Genome${view}Index"), val(pairedEnd), emit: genomeAlignmentsIndices
    tuple val(sample), val(id), path("Log.final.out"), val('txt'), val("STARstats"), val(pairedEnd), optional: true, emit: stats
    tuple val(sample), val(id), path("SJ.out.tab"), val('tsv'), val("STARjunctions"), val(pairedEnd), optional: true, emit: junctions

    script:
    type = 'bam'
    view = 'Alignments'
    def prefix = "${sample}_m${params.maxMismatches}_n${params.maxMultimaps}"

    // prepare BAM @RG tag information
    // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
    def date = ""
    def readGroupList = []
    readGroupList << ["ID", "${id}"]
    readGroupList << ["PU", "${id}"]
    readGroupList << ["SM", "${sample}"]
    if ( date ) readGroupList << ["DT", "${date}"]
    if ( params.rgPlatform ) readGroupList << ["PL", "${params.rgPlatform}"]
    if ( params.rgLibrary ) readGroupList << ["LB", "${params.rgLibrary}"]
    if ( params.rgCenterName ) readGroupList << ["CN", "${params.rgCenterName}"]
    if ( params.rgDesc ) readGroupList << ["DS", "${params.rgDesc}"]

    def fqs = reads.toString().split(" ")
    pairedEnd = (fqs.size() == 2)

    def pigzCpus = Math.min(task.cpus, 4)
    def memory = (task.memory ?: 1.GB).toBytes()
    def totalMemory = (memory * 2 / 3) as long
    def threadMemory = (totalMemory / task.cpus) as long
    def halfCpus = (task.cpus > 1 ? task.cpus / 2 : task.cpus) as int

    // Setting up params and command
    def readGroup = readGroupList.collect { it.join('=') }.join(',')
    def cmd = []
    cmd << """\
        gemtools rna-pipeline -i ${genomeDir}/genome_index.gem \\
                              -a ${annotation} \\
                              -r ${genomeDir}/transcript_index.junctions.gem \\
                              -k ${genomeDir}/transcript_index.junctions.keys \\
                              -f ${reads} \\
                              --filter-max-multi-maps ${params.maxMultimaps} \\
                              --filter-max-error-events ${params.maxMismatches} \\
                              --no-bam \\
                              -t ${task.cpus} \\
                              -q ${qualityOffset} \\
                              -n ${prefix}""".stripIndent()
    cmd << "pigz -p ${task.cpus} -dc ${prefix}.filtered.map.gz \\"
    cmd << """\
        | gem-2-sam -T ${task.cpus} \\
            -I ${genomeDir}/genome_index.gem \\
            -q offset-${qualityOffset} \\
            -l \\
            --read-group ${readGroup} \\""".stripIndent()
    if ( ! pairedEnd ) {
        cmd << "    --expect-single-end-reads \\"
        cmd << '''\
            | awk 'BEGIN{OFS=FS="\\t"}$0!~/^@/{split("1_2_8_32_64_128",a,"_");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}' \\'''.stripIndent()
    } else {
        cmd << "    --expect-paired-end-reads \\"
    }
    cmd << "| samtools view -@ ${task.cpus} -Sb - \\"
    cmd << """\
        | samtools sort -@ ${task.cpus} \\
                        -m ${threadMemory} \\
                        - \\
                        -T . \\
                        -o ${prefix}_toGenome.bam""".stripIndent()
    cmd << "samtools index ${prefix}_toGenome.bam"
    cmd.join('\n')

}