include { fastaIndex } from "../modules/fastaIndex/${params.fastaIndexTool}"
include { contig } from "../modules/contig/${params.contigTool.toLowerCase()}"
include { bigwig } from "../modules/bigwig/${params.bigwigTool.toLowerCase()}"

doBigwig = ( 'bigwig' in params.stepList )
doContig = ( 'contig' in params.stepList )

workflow signal {
    take:
      genome
      genomeAlignments

    main:
      if ( doBigwig || doContig ) {
        fastaIndex( genome )
      }
      if ( doBigwig ) {
        bigwig( fastaIndex.out, genomeAlignments )
      }
      if ( doContig ) {
        contig( fastaIndex.out, genomeAlignments )
      }

    emit:
      contigs = ( doContig ) ? contig.out : Channel.empty()
      bigwigs = ( doBigwig ) ? bigwig.out.transpose() : Channel.empty()
}