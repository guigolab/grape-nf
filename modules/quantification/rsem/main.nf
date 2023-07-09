params.rsemVersion = '1.2.21'
params.container = "grapenf/quantification:rsem-${params.rsemVersion}"
params.rsemSkipCi = false
params.rsemPlotModel = false

process index {

    tag "${genome.simpleName}-${annotation.simpleName}"
    container params.container

    input:
    path(genome)
    path(annotation)

    output:
    path('txDir')

    script:
    def cmd = []
    def genomeFile = genome.name
    def annotationFile = annotation.name

    cmd << "mkdir txDir"
    if ( genome.extension in params.comprExts ) {
        genomeFile = genome.baseName
        cmd << """\
            mkfifo ${genomeFile}
            zcat ${genome} > ${genomeFile} &
        """.stripIndent()
    }
    if ( annotation.extension in params.comprExts ) {
        annotationFile = annotation.baseName
        cmd << """\
            mkfifo ${annotationFile}
            zcat ${annotation} > ${annotationFile} &
        """.stripIndent()
    }
    cmd << """\
        rsem-prepare-reference --gtf ${annotationFile} \\
                               ${genomeFile} \\
                               txDir/RSEMref""".stripIndent()
    cmd.join('\n')
}

process quantify {

    tag "${sample}"
    container params.container

    input:
    path(quantRef)
    tuple val(sample), val(id),  path(bam), val(type), val(view), val(pairedEnd), val(readStrand)

    output:
    tuple val(sample), val(id), path("*isoforms*"), val(type), val(viewTx), val(pairedEnd), val(readStrand), emit: isoforms
    tuple val(sample), val(id), path("*genes*"), val(type), val(viewGn), val(pairedEnd), val(readStrand), emit: genes

    script:
    prefix = "${sample}"
    type = 'tsv'
    viewTx = "TranscriptQuantifications"
    viewGn = "GeneQuantifications"
    def memory = (task.memory ?: 1.GB).toMega()
    def forwardProb = null

    switch (readStrand) {
        case ~/(ANTI|MATE2_)SENSE/:
            forwardProb = '0'
            break
        case ~/(MATE1_)?SENSE/:
            forwardProb = '1'
            break
    }
    
    def cmd = []
    cmd << "sambamba sort -t ${task.cpus} -m ${memory}MB -N -M -l 0 -o - ${bam} \\"
    cmd << """\
        | rsem-calculate-expression -p ${task.cpus} \\
                                    --bam \\
                                    --seed 12345 \\
                                    --estimate-rspd  \\
                                    --no-bam-output \\""".stripIndent()
    if ( pairedEnd ) {
        cmd << """\
                            --paired-end \\"""
    }
    if ( forwardProb ) {
        cmd << """\
                            --forward-prob ${forwardProb} \\"""
    }
    if ( ! params.rsemSkipCi ) {
        cmd << """\
                            --calc-ci \\
                            --ci-memory ${memory} \\"""
    }
    cmd << """\
                            - \\
                            ${quantRef}/RSEMref \\
                            ${prefix}"""
    if ( params.rsemPlotModel ) {
        cmd << "rsem-plot-model ${prefix} ${prefix}.pdf"
    }
    cmd.join('\n')
}

process makeTable {
  tag "${feature}_${unit}s"
  container 'quay.io/biocontainers/python:3.10'

  input:
    each unit
    tuple val(feature), path(results, stageAs: 'results/*')
  
  output:
    path( "${feature}_${unit}s_table.tsv")

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
                csvfile = csv.DictReader(fd, delimiter='\\t')
                for line in csvfile:
                    element_id = line[feature_id]
                    element_value_list = []
                    v = line["${unit}"]
                    element_value = float(v) if v != 'NA' else 'NA'
                    element_value_list += [str(element_value)]
                    d.setdefault(element_id, {}).setdefault(sample, ','.join(element_value_list))

        f = open('${feature}_${unit}s_table.tsv','w')
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