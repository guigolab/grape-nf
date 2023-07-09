include { quantify; makeTable } from "${baseDir}/modules/quantification/flux"

params.fluxTableUnits = [ 'RPKM', 'reads' ]
def doQuantify = ( 'quantification' in params.stepList ) 

workflow quantification {
  take:
    genome
    annotation
    genomeAlignments
    transcriptomeAlignments
  main:
    if ( doQuantify ) {
      quantify( annotation,  genomeAlignments )
      geneResults = quantify.out.genes
        .map { it[2] }
        .collect()
        .map { [ 'genes', it ] }
      txResults = quantify.out.isoforms
        .map { it[2] }
        .collect()
        .map { [ 'isoforms', it ] }
      makeTable(params.fluxTableUnits, geneResults.mix(txResults))
    }
  emit:
    genes = ( doQuantify ) ? quantify.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  quantify.out.isoforms : Channel.empty()
    tables = ( doQuantify ) ? makeTable.out : Channel.empty()
}