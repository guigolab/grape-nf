include { index;  validateBam; calculateExpression } from "../../modules/quantification/rsem"

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
        validateBam ( transcriptomeAlignments )
        calculateExpression(
            quantificationIndex,
            transcriptomeAlignments.join(validateBam.out).map {
                it[0..-2] + [ it[-1].contains("not valid") ]
            }
        )
    }
  emit:
    genes = ( doQuantify ) ? calculateExpression.out.genes : Channel.empty()
    isoforms = ( doQuantify ) ?  calculateExpression.out.isoforms : Channel.empty()
}
