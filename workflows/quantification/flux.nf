include { quantify } from "../../modules/quantification/flux"

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
    }
  emit:
    genes = ( doQuantify ) ? quantify.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  quantify.out.isoforms : Channel.empty()
}
