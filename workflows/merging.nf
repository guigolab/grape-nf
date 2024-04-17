include { markdup } from "../modules/markdup/${params.markdupTool}"
include { sortBam as sortTranscriptome } from "../modules/sortBam/${params.mergeBamTool}"
include { mergeBam as mergeGenome } from "../modules/mergeBam/${params.mergeBamTool}"
include { mergeBam as mergeTranscriptome } from "../modules/mergeBam/${params.mergeBamTool}"

def markDuplicates = params.markDuplicates || params.removeDuplicates

workflow merging {
    take:
      genomeAlignments
      transcriptomeAlignments

    main:

      genomeAlignments
        .groupTuple(by: [0, 3, 4, 5])
        .branch {
          merge: it[2].size() > 1
          single: it[2].size() == 1
        }.set { genomeAlignments }

      mergeGenome( genomeAlignments.merge )

      genomeAlignmentsOutput = mergeGenome.out.mix( genomeAlignments.single.transpose() )
      
      if ( markDuplicates ) {
        markdup( genomeAlignmentsOutput )
        genomeAlignmentsOutput = markdup.out.dedupedAlignments
      }

      transcriptomeAlignments
        .groupTuple(by: [0, 3, 4, 5])
        .branch {
          merge: it[2].size() > 1
          single: it[2].size() == 1
        }.set { transcriptomeAlignments }
      
      sortTranscriptome( transcriptomeAlignments.merge.transpose() )
      mergeTranscriptome( sortTranscriptome.out[0].groupTuple(by: [0, 3, 4, 5]) )

    emit:
      genomeAlignments = genomeAlignmentsOutput
      transcriptomeAlignments = mergeTranscriptome.out.mix( transcriptomeAlignments.single.transpose() )
}