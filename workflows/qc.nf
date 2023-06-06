include { parseJSON } from '../modules/functions'
include { inferExp } from '../modules/inferExp'
include { bamStats } from '../modules/bamStats'

workflow QC {
  take:
    genomeAlignments
    transcriptomeAlignments
  main:
    inferExp(params.annotation, genomeAlignments)
    bamStats(params.annotation, genomeAlignments)

    inferExp.out.map {
      d = parseJSON(it[-1])
      it[0..1] + [ d.paired, d.exp ]
    }. set { inferExpOut }

    genomeAlignments
      .map {
        it[0..-2]
      }
      .join( inferExpOut, by: [0,1] )
      .set { genomeAlignmentsInferred }

    transcriptomeAlignments
      .map {
        it[0..-2]
      }
      .join( inferExpOut, by: [0,1] )
      .set { transcriptomeAlignmentsInferred }

    bamStats.out
      .map {
        it[0..-2]
      }
      .join( inferExpOut, by: [0,1] )
      .set { bamStatsInferred }
  emit:
    genomeAlignments = genomeAlignmentsInferred
    transcriptomeAlignments = transcriptomeAlignmentsInferred
    bamStats = bamStatsInferred
}