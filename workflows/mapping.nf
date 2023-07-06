include { index; map } from "../modules/mapping/${params.mappingTool.toLowerCase()}"
workflow mapping {
  take:
    genome
    annotation
    input
  main:
    if ( ! params.genomeIndex ) {
      genomeIndex = index( genome, annotation )
    } else {
      genomeIndex = file(params.genomeIndex)
      log.info "[INFO] Using precomputed genome index ${genomeIndex}"
      log.info ""
    }
    map( annotation, genomeIndex, input )
  emit:
    genomeAlignments = map.out.genomeAlignments
    transcriptomeAlignments = map.out.transcriptomeAlignments
    genomeAlignmentsIndices = map.out.genomeAlignmentsIndices
    stats = map.out.stats
    junctions = map.out.junctions
}