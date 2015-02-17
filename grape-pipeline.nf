#!/bin/env nextflow 
/*
 * Copyright (c) 2014, Centre for Genomic Regulation (CRG)
 * Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
 *
 * This file is part of the Blueprint RNAseq pipeline.
 *
 * The Blueprint RNAseq pipeline is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Set default values for params

params.steps = 'mapping,bigwig,contig,quantification'
//params.tmpDir = (System.env.TMPDIR != null ? true : false)
params.maxMismatches = 4
params.maxMultimaps = 10
params.bamStats = false
//params.countElements = [] 
//params.fluxMem = '3G'
//params.fluxProfile = false 
if (params.chunkSize) params.chunkSize = params.chunkSize as int
if (params.sjOverHang) params.sjOverHang = params.sjOverHang as int else params.sjOverHang = 100
if (!(params.wigRefPrefix)) params.wigRefPrefix = 'chr'

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

//print usage
if (params.help) {
    log.info ''
    log.info 'G R A P E ~ RNA-seq Pipeline'
    log.info '----------------------------'
    log.info 'Run the GRAPE RNA-seq pipeline on a set of data.'
    log.info ''
    log.info 'Usage: '
    log.info '    grape-pipeline.nf -i INDEX_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...'
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --index INDEX_FILE                  Index file.'
    log.info '    --genome GENOME_FILE                Reference genome file(s).'
    log.info '    --annotation ANNOTAION_FILE         Reference gene annotation file(s).'
    log.info '    --steps STEP[,STEP]...              The steps to be executed within the pipeline run. Possible values: "mapping", "bigwig", "contig", "quantification". Default: all'
//    log.info '    --tmp-dir                           Specify the temporary folder to be used as a scratch area.'
//    log.info '                                        Default: "$TMPDIR" if the environment variable is defined, "-" otherwise.'
//    log.info '    --chunk-size CHUNK_SIZE             The number of records to be put in each chunk when splitting the input. Default: no split'
//    log.info '    --paired-end                        Specify whether the data is paired-end. Default: "auto".'
    log.info '    --error-strategy ERROR_STRATEGY     Specify how an error condition is managed by the pipeline processes. Possible values: ignore, retry'
    log.info '                                        Default: the entire pipeline  terminates if a process returns an error status.'
//    log.info '    --max-read-length READ_LENGTH       The maximum read length (used to compute the transcriptomes). Default: "auto".'
    log.info '    --max-mismatches THRESHOLD          Set maps with more than THRESHOLD error events to unmapped. Default "4".'
    log.info '    --max-multimaps THRESHOLD           Set multi-maps with more than THRESHOLD mappings to unmapped. Default "10".'
//    log.info '    --filter-intron-length THRESHOLD    Filter multimaps preferring ones with intron length > THRESHOLD'
//    log.info '    --filter-block-length THRESHOLD     Filter multimaps preferring ones with block length > THRESHOLD'
//    log.info '    --filter-level LEVEL                Reduce multimaps using the specified uniqueness level.'
//    log.info '    --read-strand READ_STRAND           Directionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default "auto".'
    log.info '    --rg-platform PLATFORM              Platform/technology used to produce the reads for the BAM @RG tag.'
    log.info '    --rg-library LIBRARY                Sequencing library name for the BAM @RG tag.'
    log.info '    --rg-center-name CENTER_NAME        Name of sequencing center that produced the reads for the BAM @RG tag.'
    log.info '    --rg-desc DESCRIPTION               Description for the BAM @RG tag.'
//    log.info '    --flux-mem MEMORY                   Specify the amount of ram the Flux Capacitor can use. Default: "3G".'
//    log.info '    --flux-profile                      Specify whether the Flux Capacitor profile file should be written. Default: "false".'
//    log.info '    --count-elements ELEMENTS           A comma separated list of elements to be counted by the Flux Capacitor.'
//    log.info '                                        Possible values: INTRONS, SPLICE_JUNCTIONS. Default: "none".'
//    log.info '    --loglevel LOGLEVEL                 Log level (error, warn, info, debug). Default "info".'
    exit 1
}

if (!params.genome) {
    exit 1, "Genome file not specified"
}

if (!params.annotation) {
    exit 1, "Annotation file not specified"
}

log.info ""
log.info "B L U E P R I N T ~ RNA Pipeline"
log.info ""
log.info "General parameters"
log.info "------------------"
log.info "Index file                      : ${params.index}"
log.info "Genome                          : ${params.genome}"
log.info "Annotation                      : ${params.annotation}"
log.info "Pipeline steps                  : ${pipelineSteps.join(" ")}"
//log.info "Input chunk size                : ${params.chunkSize != null ? params.chunkSize : 'no split'}"
log.info "Error strategy                  : ${params.errorStrategy != null ? params.errorStrategy : 'default'}"
//log.info "Use temporary folder      : ${params.tmpDir}"
log.info ""
if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Max mismatches                  : ${params.maxMismatches}"
    log.info "Max multimaps                   : ${params.maxMultimaps}"
//    log.info "Max read length                 : ${params.maxReadLength != null ? params.maxReadLength : 'auto'}"
//    log.info "Read strandedness               : ${params.readStrand != null ? params.readStrand : 'auto'}"
//    log.info "Paired-end                      : ${params.pairedEnd != null ? params.pairedEnd : 'auto'}"
    log.info "Produce BAM stats               : ${params.bamStats}"
    if ( params.rgPlatform ) log.info "Sequencing platform             : ${params.rgPlatform}"  
    if ( params.rgLibrary ) log.info "Sequencing library              : ${params.rgLibrary}"  
    if ( params.rgCenterName ) log.info "Sequencing center               : ${params.rgCenterName}" 
    if ( params.rgDesc ) log.info "@RG Descritpiton                : ${params.rgDesc}" 
    log.info ""
}
if ('bigwig' in pipelineSteps) {
    log.info "Bigwig parameters"
    log.info "-----------------"
    log.info "References prefix               : ${params.wigRefPrefix != null ? params.wigRefPrefix : 'all'}"
    log.info ""
}
//if ('quantification' in pipelineSteps) {
//    log.info "Quantification parameters"
//    log.info "-------------------------"
//    log.info "Additional quantified elements  : ${params.countElements.size()==0 ? 'NONE' : params.countElements.join(" ")}"
//    log.info "Memory                          : ${params.fluxMem}"
//    log.info "Create profile file             : ${params.fluxProfile}"
//    log.info ""
//}

genomes=params.genome.split(',').collect { file(it) }
annos=params.annotation.split(',').collect { file(it) }

index = params.index ? file(params.index) : System.in
input_files = Channel.create()
input_chunks = Channel.create()

data = ['samples': [], 'ids': []]
merge = false
input = Channel
    .from(index.readLines())
    .map {
        line -> [ line.split()[0], line.split()[1], file(line.split()[2]), line.split()[3], line.split()[4] ]
    }

if (params.chunkSize) {
    merge = true
    input = input.splitFastq(by: params.chunkSize, file: true, elem: 2)
}

input.subscribe onNext: { 
        sample, id, path, type, view ->
        items = "sequencing runs"
        if( params.chunkSize ) {
            items = "chunks         "
            id = id+path.baseName.find(/\..+$/)
        } 
        input_chunks << tuple(sample, id, path, type, view)
        data['samples'] << sample
        data['ids'] << id }, 
    onComplete: { 
        ids=data['ids'].unique().size()
        samples=data['samples'].unique().size()
        if (ids != samples) merge=true
        log.info "Dataset information"
        log.info "-------------------"
        log.info "Number of sequenced samples     : ${samples}"
        log.info "Number of ${items}       : ${ids}" 
        log.info "Merging                         : ${ merge ? 'by sample' : 'none' }"
        log.info ""
        input_chunks << Channel.STOP 
    }

input_bams = Channel.create()
bam = Channel.create()

input_chunks
    .groupBy {
        sample, id, path, type, view -> id 
    }
    .flatMap ()
    .choice(input_files, input_bams) { it -> if ( it.value[0][3] == 'fastq' ) 0 else if ( it.value[0][3] == 'bam' ) 1}
    
input_files = input_files.map {        
    [it.key, it.value[0][0], it.value.collect { sample, id, path, type, view -> path }, fastq(it.value[0][2]).qualityScore()]
}

// id, sample, type, view, "${id}${prefix}.bam", pairedEnd
input_bams.map {
    [it.key, it.value[0][0], it.value[0][3], it.value[0][4], it.value.collect { sample, id, path, type, view -> path }, true].flatten()
}
.subscribe onNext: {
    bam << it
}, onComplete: {}

Genomes = Channel.create()
Annotations = Channel.create()
Channel.from(genomes)
    .merge(Channel.from(annos)) {
        g, a -> [g.name.split("\\.", 2)[0], g, a]
    }
    .separate(Genomes, Annotations) {
        a -> [[a[0],a[1]], [a[0],a[2]]]
    }

(Genomes1, Genomes2, Genomes3) = Genomes.into(3)
(Annotations1, Annotations2, Annotations3, Annotations4, Annotations5, Annotations6, Annotations7) = Annotations.into(7)

pref = "_m${params.maxMismatches}_n${params.maxMultimaps}"

if ('contig' in pipelineSteps || 'bigwig' in pipelineSteps) {
    process fastaIndex {
        input:
        set species, file(genome) from Genomes1
        set species, file(annotation) from Annotations1
    
        output:
        set species, file("${genome}.fai") into FaiIdx
    
        script:
        def command = ""
        
        command += "samtools faidx ${genome}"
        
        return command
    }
} else {
    FaiIdx = Channel.just(Channel.STOP)
}

(FaiIdx1, FaiIdx2) = FaiIdx.into(2)
 
if ('mapping' in pipelineSteps) {
    process index {
    
        input:
        set species, file(genome) from Genomes2
        set species, file(annotation) from Annotations7
    
        output:
        set species, file("genomeDir") into GenomeIdx
    
        script:
        def command = ""
       
        command += "mkdir genomeDir\n"
        command += "STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles ${genome} --sjdbGTFfile ${annotation} --sjdbOverhang ${params.sjOverHang}" 
        
        return command
    }
     
    (GenomeIdx1, GenomeIdx2) = GenomeIdx.into(2)
    
    process t_index {
    
        input:
        set species, file(genome) from Genomes3
        set species, file(annotation) from Annotations2
    
        output:
        set species, file('txDir') into TranscriptIdx
    
        script:
        def command = ""
   
        command += "mkdir txDir\n"
        command += "rsem-prepare-reference --no-polyA --gtf ${annotation} ${genome} txDir/RSEMref"
    
        return command
    }

    process mapping {
    
        input:
        set id, sample, file(reads), qualityOffset from input_files
        set species, file(annotation) from Annotations3.first()
        set species, file(genomeDir) from GenomeIdx2.first()
    
        output:
        set id, sample, type, view, file("${id}${prefix}.bam"), pairedEnd into bam
    
        script:
        type = 'bam'
        view = 'Alignments'
        prefix = pref
        
        // prepare BAM @RG tag information
        // def date = new Date().format("yyyy-MM-dd'T'HH:mmZ", TimeZone.getTimeZone("UTC"))
        def date = ""
        def readGroup = []
        readGroup << "ID:${id}" 
        readGroup << "PU:${id}" 
        readGroup << "SM:${sample}" 
        if ( date ) readGroup << "DT:${date}"
        if ( params.rgPlatform ) readGroup << "PL:${params.rgPlatform}"
        if ( params.rgLibrary ) readGroup << "LB:${params.rgLibrary}"
        if ( params.rgCenterName ) readGroup << "CN:${params.rgCenterName}"
        if ( params.rgDesc ) readGroup << "DS:${params.rgDesc}"

        def command = ""
        
        fqs = reads.toString().split(" ")
        pairedEnd = false
        if (fqs.size() == 2) pairedEnd = true 
   
        command += "STAR --runThreadN ${task.cpus} --genomeDir ${genomeDir} --readFilesIn ${reads} --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD"
        command += " --outFilterMultimapNmax ${params.maxMultimaps}   --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.0${params.maxMismatches}"
        command += " --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --readFilesCommand pigz -p${task.cpus} -dc"
        command += " --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outSAMstrandField intronMotif"
        command += " --outSAMattrRGline ${readGroup.join(' ')}"
        command += " && mv Aligned.sortedByCoord.out.bam ${id}${prefix}.bam"
        command += " && mv Aligned.toTranscriptome.out.bam ${id}${prefix}.toTranscriptome.bam"
        command += " && samtools index ${id}${prefix}.bam"
        
        return command
    }

}

if (merge) {

    singleBam = Channel.create()
    groupedBam = Channel.create()
   
    bam.groupBy() {
        id, sample, type, view, path, pairedEnd -> sample
    }
    .flatMap()
    .map {
        [it.value.size() > 1 ? it.key : it.value[0][0], it.value[0][2], it.value[0][3], it.value.collect { id, sample, type, view, path, pairedEnd -> path }, it.value[0][5]]
    }.choice(singleBam, groupedBam) {
        it[3].size() > 1 ? 1 : 0
    }
    
    process mergeBam {
        
        input:
        set id, type, view, file(bam), pairedEnd from groupedBam
    
        output:
        set id, id, type, view, "${id}${prefix}.bam", pairedEnd into mergedBam
    
        script:
        def command = ""
        prefix = pref
    
        command += "(samtools view -H ${bam} | grep -v \"@RG\";for f in ${bam};do samtools view -H \$f | grep \"@RG\";done) > header.txt\n"
        command += "samtools merge -@ ${task.cpus} -h header.txt ${id}${prefix}.bam ${bam}"
    }
    
    bam = singleBam
    .mix(mergedBam)
    .map { 
        it.flatten() 
    }
}

(bam1, bam2) = bam
.map {
    [ it[0], it[2], it[3], it[4], it[5] ]
}
.into(2)

process inferExp {
    input:
    set id, type, view, file(bam), pairedEnd from bam1
    set species, file(annotation) from Annotations4.first()

    output:
    set id, stdout into bamStrand

    script:
    def command = ""
    prefix = pref
    def genePred = "${annotation.name.split('\\.', 2)[0]}.genePred"
    def bed12 = "${annotation.name.split('\\.', 2)[0]}.bed"

    command += "set -o pipefail\n"
    command += "gtfToGenePred ${annotation} -allErrors -ignoreGroupsWithoutExons ${genePred} 2> ${genePred}.err\n"
    command += "genePredToBed ${genePred} ${bed12}\n" 
    command += "${baseDir}/bin/infer_experiment.py -i ${bam} -r ${bed12} 2> infer_experiment.log | tr -d '\\n'"
}


allBams = bam2.mix(bamStrand)
.groupBy()
.flatMap {
    bam -> bam.collect { [ it.value[0], it.value[1][-1] ].flatten() } 
}

(bam1, bam2, bam3, out) = allBams.into(4)

if (!('bigwig' in pipelineSteps)) bam1 = Channel.just(Channel.STOP)
if (!('contig' in pipelineSteps)) bam2 = Channel.just(Channel.STOP)
if (!('quantification' in pipelineSteps)) bam3 = Channel.just(Channel.STOP)

process bigwig {
    
    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam1
    set species, file(genomefai) from FaiIdx1.first()
    
    output:
    set id, type, views, file('*.bw') into bigwig

    script:
    views = []
    view = 'Signal'
    type = "bigwig"
    strand = ['unstranded': '.raw']
    bedGraphs = ['.Unique', '.UniqueMultiple']
    
    def command = ''
    command += "mkdir Signal\n"
    command += "STAR --runThreadN ${task.cpus} --runMode inputAlignmentsFromBAM --inputBAMfile ${bam} --outWigType bedGraph"
    if (readStrand != 'NONE') {
        strand = ['strand+': '.plusRaw','strand-': '.minusRaw'] 
        command += " --outWigStrand Stranded"
    } else {
        command += " --outWigStrand Unstranded"
    }
    command += " --outFileNamePrefix ./Signal/ --outWigReferencesPrefix ${params.wigRefPrefix}\n"
    
    bedGraphs.each( { bg ->
        strand.eachWithIndex( { str, istr ->
            istr+=1
            command += "bedGraphToBigWig Signal/Signal${bg}.str${istr}.out.bg"
            if ( params.wigRefPrefix ) {
                command += " <(grep -P '^${params.wigRefPrefix}' ${genomefai})"
            } else {
                command += " ${genomefai}"
            }
            command += " ${id}${bg}${str.value}.bw\n"
            views << "${str.value[1..-1].capitalize()}${view}"
        } )
    } )

    return command

}

bigwig = bigwig.reduce([:]) { files, tuple ->
    (id, type, view, path) = tuple     
    if (!files) files = []
    paths = path.toString().replaceAll(/[\[\],]/,"").split(" ")
    (1..paths.size()).each { files << [id, type, view[it-1], paths[it-1]] }
    return files
}
.flatMap()

process contig {

    input:
    set id, type, view, file(bam), pairedEnd, readStrand from bam2
    set species, file(genomefai) from FaiIdx2.first()

    output:
    set id, type, view, file('*_contigs.bed') into contig

    script:
    type = 'bed'
    view = 'Contigs'
    def command = "" 
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (readStrand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        if (pairedEnd) mateBit = (readStrand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${task.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${task.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam ${bam}\n"
    }

    command += "bamtools filter -tag NH:1 -in ${bam} -out tmp.bam\n"
    command += "mv -f tmp.bam ${bam}\n"

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${id}${it.value}.bedgraph\n"
    } )

    if (strand.size() == 2) {
        command += "contigsNew.py --chrFile ${genomefai}"
        strand.each( {
            command += " --file${it.value.substring(1,2).toUpperCase()} ${id}${it.value}.bedgraph"
        } )
        command += " | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\\t\")} print s}'"
        command += " > ${id}_contigs.bed"
    } else {
        command += "bamToBed -i ${bam} | sort -k1,1 -k2,2n"
        command += " | mergeBed"
        command += " > ${id}_contigs.bed"
    }

    return command

}

process quantification {

    input:    
    set id, type, view, file(bam), file(bai), pairedEnd, readStrand from bam3.map { [it[0],it[1],it[2],file("${it[3].toAbsolutePath().toString().replace('.bam','.toTranscriptome.bam')}"),file("${it[3].toAbsolutePath()}.bai"),it[4],it[5]] }
    set species, file(txDir) from TranscriptIdx.first()

    output:
    set id, type, viewTx, file("Quant.genes.results") into isoforms
    set id, type, viewGn, file("Quant.isoforms.results") into genes

    script:
    type = "gtf"
    viewTx = "Transcript${txDir.name.replace('.gtf','').capitalize()}"
    viewGn = "Gene${txDir.name.replace('.gtf','').capitalize()}"
    prefix = pref
    def command = ""

    command += "cat <( samtools view -H ${bam} )  <( samtools view -@ ${task.cpus} ${bam}"
    if (pairedEnd) command += " | paste -d ' ' - -"
    command += " | sort -S ${(task.memory.toBytes()/task.cpus) as long} -T ./"
    if (pairedEnd) command += " | tr ' ' '\\n'"
    command += " ) | samtools view -@ ${task.cpus} -bS - > tmp.bam"
    command += "&& mv tmp.bam ${bam}\n"

    command += "rsem-calculate-expression --bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"
    command += " -p ${task.cpus} --ci-memory ${task.memory.toMega()}" 

    if (pairedEnd) command += " --paired-end"
    if (readStrand != "NONE") command += " --forward-prob 0"

    command += " ${bam} ${txDir}/RSEMref Quant"
    command += " && rsem-plot-model Quant Quant.pdf"

    return command
}

out.mix(bigwig, contig, isoforms, genes).collectFile(name: "pipeline.db", newLine: true) {
    [it[3], it[0], it[1], it[2]].join("\t")
}
.subscribe {
    def msg = "Output files db -> ${it}"
    log.info ""
    log.info "-" * msg.size()
    log.info "Pipeline run completed."
    log.info ""
    log.info msg
    log.info "-" * msg.size()
}
