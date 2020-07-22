#!/bin/bash

CHR_RE=$'chr([1-3]|[XY])$'
CHR_READS=450
SKIP_CHR_READS=0
READS=2000
BAM_STRIP="_m4_n10_toGenome"
OUTPUT_DIR=$PWD

get_genome() {
    local genome=$1
    local prefix=$2
    paste - - <$genome| awk '$1~/^>'$CHR_RE'/{print $1,$2}' OFS="\n" > $prefix.fa
}

get_annotation() {
    local anno=$1
    local prefix=$2
    awk '$1~/^'$CHR_RE'/' $anno > $prefix.gtf
}

get_reads() {
    local bam=$1
    local reads=$2

    # AWK script
    read -r -d '' AWK_STR <<-'EOF'
    $3 ~ re {
        if(a[$3] == chrReads) {
            next
        }
        if (!reads[$1]) {
            reads[$1]++
            a[$3]++;
            tot++;
            print $1
        }
    }
    $3 == "*" && tot < totalReads {
        if (!reads[$1]) {
            reads[$1]++
            tot++
            print $1
        }
    }
EOF
    sambamba view -F 'not secondary_alignment' $f | awk -v re="$CHR_RE" -v chrReads="$CHR_READS" -v totalReads="$READS" "$AWK_STR" | sort > "$reads"
}

write_fastq() {
    local reads=$1
    local prefix=$2

    # AWK script
    read -r -d '' AWK_STR <<-'EOF'
    func record(i) {
        return "@"$1" "i":N:0:" OFS $10 OFS "+" OFS $11
    }
    NR == FNR  {
        reads[$1]++
    }
    NR > FNR {
        if ($0 ~ /^@/) {
            print
        }
        if ($1 in reads) {
            print
        }
    }
EOF
    samtools sort -n $f \
    | sambamba view -h -F 'not secondary_alignment' /dev/stdin \
    | awk -F"\t" "$AWK_STR" OFS="\t" "$reads" - \
    | samtools fastq -N -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -
}

# Parse cli options
POSITIONAL=()
while [[ $# -gt 0 ]]
do
case $1 in
    --re)
    CHR_RE=$"$2"
    shift
    shift
    ;;
    --skip)
    SKIP_CHR_READS=$2
    shift
    shift
    ;;
    --output-dir|-C)
    OUTPUT_DIR=$2
    shift
    shift
    ;;
    *)
    POSITIONAL+=("$1")
    shift
    ;;
esac
done
set -- ${POSITIONAL[@]}

# Prepare input data
tmpdir=$(mktemp -d)
trap "rm -r $tmpdir" EXIT
echo "Temp dir: $tmpdir" >&2

while read f
do
    outname=test$((++i))
    name=$(basename $f $BAM_STRIP.bam)
    echo "Reading $name" >&2
    get_reads $f $tmpdir/$name.reads
    write_fastq $tmpdir/$name.reads $OUTPUT_DIR/$outname
done
