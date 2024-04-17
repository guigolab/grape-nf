include { index; map } from "../modules/mapping/${params.mappingTool.toLowerCase()}"
include { sortBam } from "../modules/sortBam/${params.sortBamTool.toLowerCase()}"
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
    sortBam( map.out.genomeAlignments )
  emit:
    genomeAlignments = sortBam.out[0]
    transcriptomeAlignments = map.out.transcriptomeAlignments
    genomeAlignmentsIndices = sortBam.out[1]
    stats = map.out.stats
    junctions = map.out.junctions
}