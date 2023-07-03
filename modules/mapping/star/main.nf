params.starVersion = "2.4.0j"
params.container = "grapenf/mapping:star-${params.starVersion}"
params.mappingSortTool = 'sambamba'

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
        cmd << "pigz -p ${pigzCpus} -dc ${genome} > ${genomeFile}"
    }
    if ( compressedAnnotation ) {
        annotationFile = annotation.baseName
        cmd << "pigz -p ${pigzCpus} -dc ${annotation} > ${annotationFile}"
    }
    cmd << """\
        STAR --runThreadN ${task.cpus} \\
             --runMode genomeGenerate \\
             --limitGenomeGenerateRAM ${memory} \\
             --genomeDir genomeDir \\
             --genomeFastaFiles ${genomeFile} \\
             --sjdbGTFfile ${annotationFile} \\
             --sjdbOverhang ${params.sjOverHang}""".stripIndent()
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

    def txQuant = ( params.quantificationTool.toLowerCase() == 'rsem' )
    def pigzCpus = Math.min(task.cpus, 4)
    def memory = (task.memory ?: 1.GB).toBytes()
    def totalMemory = (memory * 2 / 3) as long
    def threadMemory = (totalMemory / task.cpus) as long
    def halfCpus = (task.cpus > 1 ? task.cpus / 2 : task.cpus) as int

    // Setting up params and command
    def readGroup = readGroupList.collect { it.join(':') }.join(' ')
    def maxMismatchesParam = 'outFilterMismatchNoverReadLmax'
    if ( params.starVersion.startsWith("2.3") ) {
        maxMismatchesParam = 'outFilterMismatchNoverLmax'
    }
    def outSAMtype = ['BAM', 'SortedByCoordinate']
    def sortParams = [ 
        "     --limitBAMsortRAM ${totalMemory}",
        "mv Aligned.sortedByCoord.out.bam ${prefix}_toGenome.bam"
    ]
    if ( params.mappingSortTool in params.mappingSortTools ) {
        outSAMtype[1] = 'Unsorted'
        switch(params.mappingSortTool) {
            case 'samtools':
                sortParams = [
                    "     --outStd BAM_Unsorted \\",
                    """\
                    | samtools sort -@ ${halfCpus} \\
                        -m ${threadMemory} \\
                        - \\
                        ${prefix}_toGenome""".stripIndent() ]
                break
            case 'sambamba':
                sortParams = [ 
                    "     --outStd BAM_Unsorted \\",
                    """\
                    | sambamba sort -t ${halfCpus} \\
                        -m ${totalMemory} \\
                        -o ${prefix}_toGenome.bam \\
                        /dev/stdin""".stripIndent() ]
                break
        }
    }
    def cmd = []
    cmd << """\
        STAR --runThreadN ${task.cpus} \\
             --genomeDir ${genomeDir} \\
             --readFilesIn ${reads} \\
             --outSAMunmapped Within \\
             --outFilterType BySJout \\
             --outSAMattributes NH HI AS NM MD \\
             --outFilterMultimapNmax ${params.maxMultimaps} \\
             --outFilterMismatchNmax 999 \\
             --${maxMismatchesParam} 0.0${params.maxMismatches} \\
             --alignIntronMin 20 \\
             --alignIntronMax 1000000 \\
             --alignMatesGapMax 1000000 \\
             --alignSJoverhangMin 8 \\
             --alignSJDBoverhangMin 1 \\
             --readFilesCommand pigz -p${pigzCpus} -dc \\
             --outSAMtype ${outSAMtype.join(' ')} \\
             --outSAMattrRGline ${readGroup} \\""".stripIndent()
    if ( txQuant ) {
        cmd << '     --quantMode TranscriptomeSAM \\'
    }
    if ( params.addXs ) {
        cmd << """\
             --outSAMstrandField intronMotif \\
             --outFilterIntronMotifs RemoveNoncanonical \\
        """.stripIndent()
    }
    cmd += sortParams
    if ( txQuant ) {
        cmd << "mv Aligned.toTranscriptome.out.bam ${prefix}_toTranscriptome.bam"
    }
    cmd << "${params.mappingSortTool} index ${prefix}_toGenome.bam"
    cmd.join('\n')
}