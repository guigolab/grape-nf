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
params.dbFile = 'pipeline.db'
params.genomeIndex = null
params.help = false
params.helpAll = false
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
include { checkParams; readTsv; resolveFile; printUsage; printLog; writeDB } from './modules/functions'

// Print pipeline usage if `--help` is passed
printUsage()

// Check mandatory options
checkParams(params)

// Print pipeline log
printLog()

include { mapping } from './workflows/mapping'
include { merging } from './workflows/merging'
include { quantification } from "./workflows/quantification/${params.quantificationTool.toLowerCase()}"
include { QC } from './workflows/qc'
include { signal } from './workflows/signal'
include { bamToFastq } from "./modules/bamToFastq/samtools"

workflow {
  def genome = file(params.genome)
  def annotation = file(params.annotation)

  // Read input
  def index = params.index ? file(params.index) : System.in
  def (merge, indexLines) = readTsv(index)
  params.merge = merge
  inputFastqs = Channel.empty()
  Channel.from(indexLines)
    .filter { it }  // get only non-empty lines
    .map { line ->
        def (sampleId, runId, fileName, format, readId) = line.split()
        fileName = resolveFile(fileName, index)
        [sampleId, runId, fileName, format, readId]
    }.set { inputFiles }

  inputFiles.branch {
	fastqs: it[3] == 'fastq'
	bams: it[3] == 'bam'
    }.set{ input }
    
  if ( 'mapping' in params.stepList ) {
    bamToFastq(
      input.bams.filter {
          it[4] == 'GenomeAlignments'
      }
    )
    bamToFastq.out
      .transpose()
      .mix( input.fastqs )
      .set { inputFastqs }
  }
  inputFastqs
    .groupTuple(by: [0,1,3], sort: true)
    .map {
      it << fastq(it[2][0]).qualityScore()
    }.set { mappingInput }

  mapping( genome, annotation, mappingInput )

  merging(mapping.out.genomeAlignments, mapping.out.transcriptomeAlignments)

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

  writeDB(params.dbFile, pipelineResults)
}

workflow.onComplete {
  log.info ""
  log.info "-----------------------"
  log.info "Pipeline run completed."
  log.info "-----------------------"
}
