#!/bin/bash
#!/usr/bin/env

FILE=""
FPPROB=0.001
KMERSIZE=31
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -f|--file)
            FILE="$2"
            shift # past argument
            shift # past value
            ;;
        -p|--fpprob)
            FPPROB="$2"
            shift
            shift
            ;;
        -k|--kmer)
            KMERSIZE="$2"
            shift
            shift
            ;;
        -n|--numkmers)
            n="$2"
            shift
            shift
            ;;
        *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

p=${FPPROB}
k=5

if [ -z "${n}" ]; then
    command -v ntcard >/dev/null 2>&1 || { echo >&2 "ntcard seems not to be installed (or in the path).  Aborting."; exit 1; }
    tfile=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`
    echo "Running ntcard to estimate number of distinct k-mers"
    ntcard -k ${KMERSIZE} -t 16 -p ${tfile}/card ${FILE}
    ofname="${tfile}/card_k${KMERSIZE}.hist"
    n=`head -n2 ${ofname} | cut -f 2 | tail -n 1`
    echo "Num distinct k-mers estimated as ${n}"
    rm -fr ${tfile}
else
    echo "n = ${n}"
fi

s="
define ceil(x) {
    auto os,xx;x=-x;os=scale;scale=0
    xx=x/1;if(xx>x).=xx--
    scale=os;return(-xx)
}

lgp_k = l($p) / $k
r = (-$k) / l(1 - e(lgp_k))
ceil(l(ceil($n * r))/l(2))
"
res=$(echo "$s" | bc -l)
echo "Set param as 2^M where M=${res}"
