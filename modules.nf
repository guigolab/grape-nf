pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

// Some configuration variables
mappingTool = "${config.process.'withName:mapping'.ext.tool} ${config.process.'withName:mapping'.ext.version}"
bigwigTool = "${config.process.'withName:bigwig'.ext.tool} ${config.process.'withName:bigwig'.ext.version}"
quantificationTool = "${config.process.'withName:quantification'.ext.tool} ${config.process.'withName:quantification'.ext.version}"
quantificationMode = "${config.process.'withName:quantification'.ext.mode}"
useContainers = config.docker?.enabled ? 'docker' : (config.singularity?.enabled ? 'singularity' : 'no')
errorStrategy = config.process.errorStrategy
executor = config.process.executor ?: 'local'
queue = config.process.queue

process fastaIndex {

    input:
    set species, file(genome)
    set species, file(annotation)

    output:
    set species, file { "${genome}.fai" } 

    script:
    template(task.ext.command)

}

process index {

    input:
    set species, file(genome)
    set species, file(annotation)

    output:
    set species, file("genomeDir")

    script:
    sjOverHang = params.sjOverHang
    readLength = params.readLength

    template(task.ext.command)

}

process mapping {

    input:
    set id, sample, file(reads), qualityOffset 
    set species, file(annotation) 
    set species, file(genomeDir) 

    output:
    set id, sample, type, view, file("*.bam"), pairedEnd into bam

    script:
    type = 'bam'
    view = 'Alignments'
    prefix = "${id}${pref}"
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
    readGroup = task.ext.readGroup

    fqs = reads.toString().split(" ")
    pairedEnd = (fqs.size() == 2)
    taskMemory = task.memory ?: 1.GB
    totalMemory = taskMemory.toBytes()
    threadMemory = taskMemory.toBytes()/(2*task.cpus)
    task.ext.sort = params.bamSort ?: task.ext.sort
    halfCpus = task.cpus / 2

    template(task.ext.command)

}

process txIndex {

    input:
    set species, file(genome)
    set species, file(annotation)

    output:
    set species, file('txDir')

    script:
    template(task.ext.command)

}

process mergeBam {

    input:
    set id, sample, type, view, file(bam), pairedEnd

    output:
    set id, sample, type, view, file("${prefix}.bam"), pairedEnd

    script:
    id = id.sort().join(':')
    prefix = "${sample}${pref}_to${view.replace('Alignments','')}"

    template(task.ext.command)

}


process inferExp {

    input:
    set id, sample, type, view, file(bam), pairedEnd
    set species, file(annotation) 

    output:
    set id, stdout

    script:
    prefix = "${annotation.name.split('\\.', 2)[0]}"

    template(task.ext.command)
}

process bigwig {

    input:
    set id, sample, type, view, file(bam), pairedEnd, readStrand 
    set species, file(genomeFai) 

    output:
    set id, sample, type, views, file('*.bw'), pairedEnd, readStrand

    script:
    type = "bigWig"
    prefix = "${sample}"
    wigRefPrefix = params.wigRefPrefix ?: ""
    views = task.ext.views

    template(task.ext.command)

}

process contig {

    input:
    set id, sample, type, view, file(bam), pairedEnd, readStrand
    set species, file(genomeFai)

    output:
    set id, sample, type, view, file('*.bed'), pairedEnd, readStrand

    script:
    type = 'bed'
    view = 'Contigs'
    prefix = "${sample}.contigs"

    template(task.ext.command)

}

process quantification {

    input:
    set id, sample, type, view, file(bam), pairedEnd, readStrand
    set species, file(quantRef) 

    output:
    set id, sample, type, viewTx, file("*isoforms*"), pairedEnd, readStrand //into isoforms
    set id, sample, type, viewGn, file("*genes*"), pairedEnd, readStrand //into genes

    script:
    prefix = "${sample}"
    refPrefix = quantRef.name.replace('.gtf','').capitalize()
    type = task.ext.fileType
    viewTx = "Transcript${refPrefix}"
    viewGn = "Gene${refPrefix}"
    memory = {task.memory ?: 1.GB}().toMega()

    template(task.ext.command)

}

/*
 * Given a string path resolve it against the index file location.
 * Params: 
 * - str: a string value represting the file pah to be resolved
 * - index: path location against which relative paths need to be resolved 
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str) 
  }
} 

def testResolveFile() {
  def index = file('/path/to/index')
  assert resolveFile('str', index) == file('/path/to/str')
  assert resolveFile('/abs/file', index) == file('/abs/file')
  assert resolveFile('s3://abs/file', index) == file('s3://abs/file')
}
