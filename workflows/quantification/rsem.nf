include { index;  quantify; makeTable } from "${baseDir}/modules/quantification/rsem"

params.rsemTableUnits = ['TPM','expected_count']
def doQuantify = ( 'quantification' in params.stepList )

workflow quantification {
  take:
    genome
    annotation
    genomeAlignments
    transcriptomeAlignments
  main:
    if ( doQuantify ) {
        quantificationIndex = index( genome, annotation )
        quantify( quantificationIndex, transcriptomeAlignments )

        geneResults = quantify.out.genes
          .map { it[2] }
          .collect()
          .map { [ 'genes', it ] }
        txResults = quantify.out.isoforms
          .map { it[2] }
          .collect()
          .map { [ 'isoforms', it ] }
        makeTable(params.rsemTableUnits, geneResults.mix(txResults))
    }
  emit:
    genes = ( doQuantify ) ? quantify.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  quantify.out.isoforms : Channel.empty()
}