params.fluxVersion = '1.6.1'
params.container = "grapenf/quantification:flux-${params.fluxVersion}"
params.fluxProfile = false

process quantify {

    tag "${sample}"
    container params.container

    input:
    path(annotation)
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    prefix = "${sample}"
    type = 'gtf'
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    def memory = (task.memory ?: 1.GB).toMega()
    def mode = pairedEnd ? 'PAIRED' : 'SINGLE'

    def cmd = []
    cmd << "samtools index ${bam}"
    if ( params.fluxProfile ) {
        cmd << """\
            flux-capacitor --profile \\
                        -i ${bam} \\ 
                        -a ${annotation} \\""".stripIndent()
        if ( readStrand in ['SENSE', 'ANTISENSE', 'MATE1_SENSE', 'MATE2_SENSE'] ) {
            cmd << """\
            -m ${mode}_STRANDED \\
            --read-strand ${readStrand} \\"""
        }
        cmd << "            --profile-file ${prefix}_profile.json"
    }
    cmd << """\
        flux-capacitor -i ${bam} \\
                    -a ${annotation} \\""".stripIndent()
    if ( readStrand in ['SENSE', 'ANTISENSE', 'MATE1_SENSE', 'MATE2_SENSE'] ) {
        cmd << """\
            -m ${mode}_STRANDED \\
            --read-strand ${readStrand} \\"""
    }
    cmd << "            -o ${prefix}.isoforms.gtf"
    cmd << """\
        TrtoGn_RPKM.sh -a ${annotation} \\
                    -i ${prefix}.isoforms.gtf \\
                    -o ${prefix}.genes.gff""".stripIndent()
    cmd.join('\n')
}

process makeTable {
  tag "${feature}_${unit.replaceAll(/s$/,'')}s"
  container 'quay.io/biocontainers/python:3.10'
  publishDir "${params.outDir}/quantificationTables", mode: 'copy'

  input:
    each unit
    tuple val(feature), path(results, stageAs: 'results/*')
  
  output:
    tuple val(feature), val("${unit.replaceAll(/s$/,'')}"), path("${feature}_${unit.replaceAll(/s$/,'')}s_table.tsv")

  script:
    def cmd = """\
        #!/usr/bin/env python
        import csv
        import os
        import sys

        featMap = {
            "genes": "gene",
            "isoforms": "transcript"
        }
        d = {}
        samples = set()

        files = os.listdir("results")
        feature_id = f'{featMap["${feature}"]}_id'

        for f in files:
            with open(f"results/{f}") as fd:
                sample = os.path.basename(f).split('.')[0]
                samples.add(sample)
                csvfile = csv.reader(fd, delimiter='\\t')
                for l in csvfile:
                    line = dict(tuple(i.strip().replace('"','').strip(',').split()) for i in filter(None, l[8].split(';')))
                    element_id = line[feature_id]
                    element_value_list = []
                    v = line["${unit}"]
                    element_value = float(v) if v != 'NA' else 'NA'
                    element_value_list += [str(element_value)]
                    d.setdefault(element_id, {}).setdefault(sample, ','.join(element_value_list))

        f = open('${feature}_${unit.replaceAll(/s$/,'')}s_table.tsv','w')
        f.write(feature_id+'\\t')
        f.write('\\t'.join(sorted(samples, key=lambda x: x.lower()))+'\\n')
        for element in sorted(d.keys()):
            values = d[element]
            f.write(element+'\t')
            missing_value = 'NA'
            f.write('\\t'.join(str(values.get(s, missing_value)) for s in sorted(samples, key=lambda x: x.lower()))+'\\n')
        f.close()\
    """.stripIndent()

    cmd
}