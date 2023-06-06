process fetch {
    tag { outPath.name }
    storeDir { outPath.parent }

    input:
    tuple val(sample), val(id), path(path), val(type), val(view)

    output:
    tuple val(sample), val(id), path("${outPath.name}"), val(type), val(view)

    script:
    def paths = path.split(',')
    outPath = workflow.launchDir.resolve(paths[-1])
    urls = paths.size() > 1 ? paths[0..-2].join(' ') : ''
    """
    for url in $urls; do
        if wget \${url}; then
            exit 0
        fi
    done
    """
}