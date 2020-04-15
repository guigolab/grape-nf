#!/bin/bash

CHR_RE=$'chr([1-3]|[XY])$'
CHR_READS=450
SKIP_CHR_READS=0
READS=2000
BAM_STRIP="_m4_n10_toGenome"
OUTPUT_DIR=$PWD

sambamba() {
    docker run --rm -v $PWD:$PWD -w $PWD grapenf/sambamba:0.7.1 sambamba -q "$@"
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

    local mate1="$(dirname $reads)/_pipe_${name}_1"
    local mate2="$(dirname $reads)/_pipe_${name}_2"
    trap "rm -f $mate1 $mate2" EXIT

    if [ ! -p "$mate1" ]; then
        mkfifo "$mate1"
    fi
    if [ ! -p "$mate2" ]; then
        mkfifo "$mate2"
    fi

    # AWK script
    read -r -d '' AWK_STR <<-'EOF'
    func record() {
        return "@"$1 OFS $10 OFS "+" OFS $11
    }
    NR == FNR  {
        reads[$1]++
    }
    NR > FNR {
        if ($1 in reads) {
            if (and($2,64)==64) {
                print record() > mate1
            }
            if (and($2,128)==128) {
                print record() > mate2
            }
        }
    }
EOF
    sambamba view -F 'not secondary_alignment' $f | awk -F"\t" -v mate1="$mate1" -v mate2="$mate2" "$AWK_STR" OFS="\t" "$reads" - &
    exec 3<"$mate1" \
        && cat <&3 | sort -k1,1 | tr '\t' '\n' | gzip -c > ${prefix}_1.fastq.gz \
        && exec 3<&- &
    exec 4<"$mate2" \
        && cat <&4 | sort -k1,1 | tr '\t' '\n' | gzip -c > ${prefix}_2.fastq.gz \
        && exec 4<&-
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
echo "Temp dir: $tmpdir"

while read f
do
    outname=test$((++i))
    name=$(basename $f $BAM_STRIP.bam)
    echo "Reading $name"
    get_reads $f $tmpdir/$name.reads
    write_fastq $tmpdir/$name.reads $OUTPUT_DIR/$outname
done