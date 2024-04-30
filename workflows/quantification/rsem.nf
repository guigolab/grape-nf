include { index;  quantify } from "../../modules/quantification/rsem"

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
    }
  emit:
    genes = ( doQuantify ) ? quantify.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  quantify.out.isoforms : Channel.empty()
}
