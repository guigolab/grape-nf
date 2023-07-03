#!/bin/env nextflow
/*
 * Copyright (c) 2015, Centre for Genomic Regulation (CRG)
 * Emilio Palumbo, Alessandra Breschi and Sarah Djebali.
 *
 * This file is part of the GRAPE RNAseq pipeline.
 *
 * The GRAPE RNAseq pipeline is a free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

nextflow.enable.dsl = 2

// Set default values for params
params.addXs = false
params.chunkSize = null
params.dbFile = 'pipeline.db'
params.genomeIndex = null
params.help = false
params.markDuplicates = false
params.removeDuplicates = false
params.maxMismatches = 4
params.maxMultimaps = 10
params.pairedEnd = false
params.readStrand = null
params.rgCenterName = null
params.rgDesc = null
params.rgLibrary = null
params.rgPlatform = null
params.sjOverHang = 100
params.steps = 'mapping,bigwig,contig,quantification'
params.stepList = params.steps.split(',').collect { it.trim() }
params.inferExpThreshold = 0.8
params.mappingSortTools = [ 'samtools', 'sambamba' ]
params.comprExts = ['gz', 'bz2', 'zip']

// Import functions
include { readTsv; resolveFile; printUsage; printLog } from './modules/functions'

// Print pipeline usage if `--help` is passed
printUsage()

// Check mandatory options
if (!params.genomeIndex && !params.genome) {
    exit 1, "Reference genome not specified"
}

if ('quantification' in params.stepList && !params.annotation) {
    exit 1, "Annotation not specified"
}

// Print pipeline log
printLog()

// Init I/O files
pdb = file(params.dbFile)
index = params.index ? file(params.index) : System.in
(merge, indexLines) = readTsv(index)
params.merge = merge

include { mapping } from './workflows/mapping'
include { merging } from './workflows/merging'
include { quantification } from "./workflows/quantification/${params.quantificationTool.toLowerCase()}"
include { QC } from './workflows/qc'
include { signal } from './workflows/signal'

workflow {
  def genome = file(params.genome)
  def annotation = file(params.annotation)

  Channel.from(indexLines)
    .filter { it }  // get only non-empty lines
    .map { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        fileName = resolveFile(fileName, index)
        [sampleId, runId, fileName, format, readId]
    }.tap {
        inputFiles
    }

  inputFiles
    .filter {
      it[3] == 'fastq'
    }
    .groupTuple(by: [0,1,3], sort: true)
    .map {
      it << fastq(it[2][0]).qualityScore()
    }.set { mappingInput }

  mapping( genome, annotation, mappingInput )

  inputFiles
    .filter {
      it[3] == 'bam'
    }
    .map { 
        it << params.pairedEnd
        it.flatten()
    }
    .mix(mapping.out.genomeAlignments)
    .mix(mapping.out.transcriptomeAlignments)
    .branch {
        genome: it[4] == 'GenomeAlignments'
        transcriptome: it[4] == 'TranscriptomeAlignments'
    }
    .set { mappings }

  merging(mappings.genome, mappings.transcriptome)

  QC(merging.out)

  signal(genome, QC.out.genomeAlignments)

  quantification(genome, annotation, QC.out.genomeAlignments, QC.out.transcriptomeAlignments)
  
  // Mix results
  QC.out.genomeAlignments.mix(
    QC.out.transcriptomeAlignments, 
    QC.out.bamStats, 
    signal.out.bigwigs, 
    signal.out.contigs, 
    quantification.out.isoforms, 
    quantification.out.genes
  )
  .set { pipelineResults }
  
  // Clear pipeline db
  pdb.write('')
  
  // Write results to pipeline db
  pipelineResults.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) {
    res = it[0..4]
    res << (it[5] ? 'Paired-End' : 'Single-End')
    res << it[6]
    res.join('\t')
  }

}

workflow.onComplete {
  log.info ""
  log.info "-----------------------"
  log.info "Pipeline run completed."
  log.info "-----------------------"
}
