include { txIndex as index } from '../modules/txIndex'
include { quantification as quantify } from '../modules/quantification'

def doQuantify = ( 'quantification' in params.stepList ) 

workflow quantification {
  take:
    genome
    annotation
    genomeAlignments
    transcriptomeAlignments
  main:
    quantificationIndex = annotation
    quantificationInput = genomeAlignments
    if ( doQuantify ) {
      if ( params.quantificationMode == 'Transcriptome' ) {
        quantificationIndex = index( genome, annotation )
        quantificationInput = transcriptomeAlignments
      }
      quantify( quantificationInput,  quantificationIndex )
    }
  emit:
    genes = ( doQuantify ) ? quantify.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  quantify.out.isoforms : Channel.empty()
}