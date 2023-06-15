include { index; map } from "../modules/mapping/${params.mappingTool.toLowerCase()}"
workflow mapping {
  take:
    genome
    annotation
    input
  main:
    if ( ! params.genomeIndex ) {
      index( genome, annotation )
    }
    map( annotation, index.out, input )
  emit:
    genomeAlignments = map.out.genomeAlignments
    transcriptomeAlignments = map.out.transcriptomeAlignments
    genomeAlignmentsIndices = map.out.genomeAlignmentsIndices
    stats = map.out.stats
    junctions = map.out.junctions
}