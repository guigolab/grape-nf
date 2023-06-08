def pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

process mapping {

    label "mapping"
    tag params.tag

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
    prefix = "${sample}${pref}"
    tagSuffix = "-${params.mappingTool}-${params.mappingToolVersion}"
    maxMultimaps = params.maxMultimaps
    maxMismatches = params.maxMismatches

    // prepare BAM @RG tag information
    // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
    date = ""
    readGroupList = []
    readGroupList << ["ID", "${id}"]
    readGroupList << ["PU", "${id}"]
    readGroupList << ["SM", "${sample}"]
    if ( date ) readGroupList << ["DT", "${date}"]
    if ( params.rgPlatform ) readGroupList << ["PL", "${params.rgPlatform}"]
    if ( params.rgLibrary ) readGroupList << ["LB", "${params.rgLibrary}"]
    if ( params.rgCenterName ) readGroupList << ["CN", "${params.rgCenterName}"]
    if ( params.rgDesc ) readGroupList << ["DS", "${params.rgDesc}"]
    (s,t) = params.mappingReadGroupSeparators
    readGroup = readGroupList.collect { it.join(s) }.join(t)

    fqs = reads.toString().split(" ")
    pairedEnd = (fqs.size() == 2)
    taskMemory = task.memory ?: 1.GB
    totalMemory = (taskMemory.toBytes()*2/3) as int
    threadMemory = (totalMemory/task.cpus) as int
    cpus = task.cpus
    halfCpus = (task.cpus > 1 ? task.cpus / 2 : task.cpus) as int

    command = "${params.mappingTool}-${params.mappingToolVersion.split("\\.")[0..1].join(".")}"
    switch(params.mappingTool) {
        case 'GEM':
            command += "-${pairedEnd ? 'Paired-End' : 'Single-End'}"
            break
        case 'STAR':
            command += (params.mappingSortTool ? "-"+params.mappingSortTool : '') + (params.quantificationMode ? "-"+params.quantificationMode : '') + (params.addXs ? "-XS" : '')
            break
    }
    template(command)

}