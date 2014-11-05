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

// get list of steps from comma-separated strings
pipelineSteps = params.steps.split(',').collect { it.trim() }

//print usage
if (params.help) {
    println '''
B L U E P R I N T ~ RNA Pipeline
--------------------------------
Run the RNAseq pipeline on one sample.

Usage: 
    ./blueprint.pipeline.sh -i FASTQ_FILE -g GENOME_FILE -a ANNOTATION_FILE [OPTION]...

Options:
    --index               index file.
    --genome              reference genome file.
    --annotation          reference gene annotation file.
    --mismatches          Max number of mismatches. Default "4".
    --hits                Max number of hits. Default "10".
    --quality-offset      The quality offset of the fastq files. Default: "33".
    --max-read-length     The maximum read length (used to compute the transcriptomes). Default: "150".
    --read-strand         directionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default "NONE".
    --loglevel            Log level (error, warn, info, debug). Default "info".
    --threads             Number of threads. Default "1".
    --paired-end          Specify whether the data is paired-end. Defalut: "false".
    --count-elements      A comma separated list of elements to be counted by the Flux Capacitor.
                          Possible values: INTRONS,SPLICE_JUNCTIONS. Defalut: "none".
    --help                Show this message and exit.
    --bam-stats           Run the RSeQC stats on the bam file. Default "false".
    --flux-mem            Specify the amount of ram the Flux Capacitor can use. Default: "3G".
    --tmp-dir             Specify local temporary folder to copy files when running on shared file systems.
                          Default: "$TMPDIR" if the environment variable is defined, "-" otherwise.
    --dry-run             Test the pipeline. Writes the command to the standard output.
    '''
    exit 1
}

// Check required parameters
if (!params.input) {
    exit 1, "Input file not specified"
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
log.info "Input                     : ${params.input}"
log.info "Genome                    : ${params.genome}"
log.info "Annotation                : ${params.annotation}"
log.info "Steps to be performed     : ${params.steps.replace(',',' ')}"
log.info "Use temporary folder      : ${params.tmp_dir}"
log.info "Number of cpus            : ${params.cpus}"
log.info ""
if ('mapping' in pipelineSteps) {
    log.info "Mapping parameters"
    log.info "------------------"
    log.info "Max mismatches            : ${params.mismatches}"
    log.info "Max multimaps             : ${params.hits}"
    log.info "Quality offset            : ${params.quality_offset}"
    log.info "Read length               : ${params.max_read_length}"
    log.info "Read strandedness         : ${params.read_strand}"
    log.info "Paired                    : ${params.paired_end}"
    log.info "Read group                : ${params.read_group}"
    log.info "Produce BAM stats         : ${params.bam_stats}"
    log.info ""
}
if ('flux' in pipelineSteps || 'quantification' in pipelineSteps) {
    log.info "Flux Capacitor parameters"
    log.info "-------------------------"
    log.info "Elements to be quantified : ${params.count_elements}"
    log.info "Memory                    : ${params.flux_mem}"
    log.info "Create profile file       : ${params.flux_profile}"
    log.info ""
}

genomes=params.genome.split(',').collect { file(it) }
annos=params.annotation.split(',').collect { file(it) }

index = file(params.index)

input_files = Channel
    .from(index.readLines())
    .map {
        line -> [ file(line.split()[0]), line.split()[1], line.split()[2], line.split()[3] ]
    }
    .groupBy {
        path, id, type, view ->  path.name.split("\\.", 2)[0][0..-3]
    }
    .flatMap ()
    .map {
       [it.key, it.value.collect { path, id, type, view -> path }].flatten()
    }
    .subscribe {
        println it
    }

return


process index {
    input:
    file genome_file

    output:
    file "genome_index.gem" into genome_index

    script:
    """
    gemtools index -i ${genome_file} -t ${params.cpus} -o genome_index.gem
    """
}

(genome_index1, genome_index2, genome_index3) = genome_index.into(3)


process t_index {
    input:
    file genome_index from genome_index1
    file annotation_file

    output:
    set file('tx_index.junctions.gem'), file('tx_index.junctions.keys') into tx_index

    script:
    """
    gemtools t-index -i ${genome_index} -a ${annotation_file} -m ${params.max_read_length} -t ${params.cpus} -o tx_index
    """
}


process mapping {
    input:
    set reads_name, file(read1), file(read2) from input_files
    file genome_index from genome_index2.first()
    set file(tx_index), file(tx_keys) from tx_index.first()

    output:
    set reads_name, view, "mapping.map.gz" into map

    script:
    view = 'gemUnfiltered'

    def command = "gemtools rna-pipeline -i ${genome_index} -r ${tx_index} -k ${tx_keys} -f ${read1}"
    command += " --no-stats --no-bam"
    if (!params.paired_end) {
        command += " --single-end"
    }
    command += " -t ${params.cpus} -q ${params.quality_offset} -n mapping"

    return command
}

process filter {
    input:
    set reads_name, view, gem_unfiltered from map

    output:
    set reads_name, view, "mapping_filtered.map.gz" into fmap

    script:
    view = "gemFiltered"

    def command = "gt.quality -i ${gem_unfiltered} -t ${params.cpus}"
    command += " | gt.filter --max-levenshtein-error ${params.mismatches} -t ${params.cpus}"
    command += " | gt.filter --max-matches ${params.hits} -t ${params.cpus}"
    command += " | pigz -p ${params.cpus} -c"
    command += " > mapping_filtered.map.gz"

    return command
}

(map1, map2) = fmap.into(2)

process gemStats {
    input:
    set reads_name, view, gem_filtered from map1

    output:
    set reads_name, view, "mapping_filtered.map.gz.stats" into stats

    script:
    view = "gemFilteredStats"

    command="gt.stats -i ${gem_filtered} -t ${params.cpus} -a"
    if (params.paired_end) {
        command += " -p"
    }
    command += " 2> mapping_filtered.map.gz.stats"
}

process gemToBam {
    input:
    set reads_name, view, gem_filtered from map2
    file genome_index from genome_index3.first()

    output:
    set reads_name, view, "mapping.bam" into bam

    script:
    view = "alignments"

    def command = "pigz -p ${params.cpus} -dc ${gem_filtered}"
    command += " | gem-2-sam -T ${params.cpus} -I ${genome_index} -q offset-${params.quality_offset} -l"
    if (params.read_group) {
       command += " --read-group ${params.read_group}"
    }
    if (params.paired_end) {
       command += " --expect-paired-end-reads"
    }
    else {
       command += " --expect-single-end-reads"
       command += " | awk 'BEGIN{OFS=FS=\"\t\"}\$0!~/^@/{split(\".toString()1_2_8_32_64_128\",a,\"_\");for(i in a){if(and(\$2,a[i])>0){\$2=xor(\$2,a[i])}}}{print}'"
    }

    command += " | samtools view -@ ${params.cpus} -Sb -"
    command += " | samtools sort -@ ${params.cpus} -m 4G - mapping"
    command += " && samtools index mapping.bam"

    return command
}

(bam1, bam2, bam3, bam4) = bam.into(4)

genomeFai = file(params.genome+'.fai')

process bigwig {
    input:
    set reads_name, in_view, file(bam) from bam1
    file genomeFai

    output:
    set reads_name, view, file('*.bw') into bigwig

    script:
    view = 'bigwig'
    def command = ''
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (params.read_strand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        mateBit = (params.read_strand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${params.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${params.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam mapping.bam\n"
    }

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${reads_name}${it.value}.bedgraph\n"
        command += "bedGraphToBigWig ${reads_name}${it.value}.bedgraph ${genomeFai} ${reads_name}${it.value}.bw\n"
    } )

    return command

}

process contig {
    input:
    set reads_name, in_view, file(bam) from bam2
    file genomeFai

    output:
    set reads_name, view, file('*_contigs.bed') into contig

    script:
    view = 'contig'
    def command = ''
    strand = ['': '']
    mateBit = 0
    awkCommand = 'BEGIN {OFS=\"\\t\"} {if (\$1!~/^@/ && and(\$2,MateBit)>0) {\$2=xor(\$2,0x10)}; print}'
    if (params.read_strand != 'NONE') {
        strand = ['+': '.plusRaw','-': '.minusRaw']
        mateBit = (params.read_strand =~ /MATE2/ ? 64 : 128)
    }

    if (mateBit > 0) {
        command += "samtools view -h -@ ${params.cpus} ${bam}"
        command += " | awk -v MateBit=${mateBit} '${awkCommand}'"
        command += " | samtools view -@ ${params.cpus} -Sb -"
        command += " > tmp.bam\n"
        command += "mv -f tmp.bam mapping.bam\n"
    }

    command += "bamflag -in ${bam} -out tmp.bam -m 3\n"
    command += "mv -f tmp.bam mapping.bam\n"

    strand.each( {
        command += "genomeCoverageBed "
        command += (it.key != '' ? "-strand ${it.key} ".toString() : ''.toString())
        command += "-split -bg -ibam ${bam} > ${reads_name}${it.value}.bedgraph\n"
    } )

    if (strand.size() == 2) {
        command += "contigsNew.py --chrFile ${genomeFai}"
        strand.each( {
            command += " --file${it.value.substring(1,2).toUpperCase()} ${reads_name}${it.value}.bedgraph"
        } )
        command += " | awk '{s=\"\"; for(i=1; i<=NF; i++){s=(s)(\$i)(\"\t\")} print s}'"
        command += " > ${reads_name}_contigs.bed"
    } else {
        command += "bamToBed -i ${bam} | sort -k1,1 -k2,2n"
        command += " | mergeBed"
        command += " > ${reads_name}_contigs.bed"
    }

    return command

}

process quantification {
    input:
    set reads_name, in_view, file(bam) from bam3
    file annotation_file

    output:
    set reads_name, view, file('flux.gtf') into flux

    script:
    view = 'transcript'
    def command = ""
    paramFile = file('params.flux')

    //Flux parameter file has to be in the same folder as the bam - Flux bug

    /* paramFile.append("# Flux Capacitor parameter file for ${reads_name}\n")
    annotationMapping = "AUTO"
    if (params.read_strand != "NONE") {
        paramFile.append("READ_STRAND ${params.read_strand}\n")
        annotationMapping="STRANDED"
        if (params.paired_end) {
            annotationMapping="PAIRED_${annotationMapping}"
        }
        else {
            annotationMapping="SINGLE_${annotationMapping}"
        }
    }
    paramFile.append("ANNOTATION_MAPPING ${annotationMapping}\n")
    paramFile.append("COUNT_ELEMENTS ${params.count_elements}\n")

    if (params.flux_profile) {
        paramFile.append("PROFILE_FILE profile.json\n")
        command += "flux-capacitor --profile -p ${paramFile} -i ${bam}  -a ${annotation_file}"
    }
    command += "flux-capacitor -p ${paramFile} -i ${bam} -a ${annotation_file} -o flux.gtf".toString() */

    // Workaround
    flux_params = ''
    annotationMapping = "AUTO"
    if (params.read_strand != "NONE") {
        flux_params += " --read-strand ${params.read_strand}"
        annotationMapping="STRANDED"
        if (params.paired_end) {
            annotationMapping="PAIRED_${annotationMapping}"
        }
        else {
            annotationMapping="SINGLE_${annotationMapping}"
        }
    }
    flux_params += " -m ${annotationMapping}"
    flux_params += " --count-elements ${params.count_elements}"

    if (params.flux_profile) {
        flux_params += " --profile-file profile.json"
        command += "flux-capacitor --profile ${flux_params} -i ${bam}  -a ${annotation_file} && "
    }
    command += "flux-capacitor ${flux_params} -i ${bam} -a ${annotation_file} -o flux.gtf"

    return command
}

bam4.mix(bigwig, contig, flux).subscribe {
    println it
}

//process store {
//
//    maxForks 1
//
//    input:
//    val hash
//    set reads_name, view, file(store_file) from bam4.mix(bigwig, contig, flux)
//
//    script:
//    """
//    idxtools add path=`readlink -f ${store_file}` id=${reads_name} view=${view} type=${store_file.name.split("\\.", 2)[1]} size=`cat ${store_file} | wc -c` md5sum=`md5sum ${store_file} | cut -d" " -f1`
//    """
//}



