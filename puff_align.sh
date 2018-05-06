#!/bin/bash
#!/usr/bin/env

start=`date +%s`

config='config.json'
if [[ $# -eq 1 ]]; then
	config=$1
fi

pufferfish=$(jq -r '.pufferfish' ${config})/build/src/pufferfish
salmon=$(jq -r '.salmon' ${config})/build/src/salmon
puffidx=$(jq -r '.puffidx' ${config})
fasta=$(jq -r '.fasta' ${config})
ratio1=$(jq -r '.ratio' ${config})
read1=$(jq -r '.read1' ${config})
read2=$(jq -r '.read2' ${config})
priorr=$(jq -r '.prior' ${config})
outdir=$(jq -r '.outDir' ${config})
outfile=$(jq -r '.outFile' ${config})
outfile=${outfile}_sr${ratio1}
truthfile=$(jq -r '.truthFile' ${config})
taxatruthfile=$(jq -r '.taxaTruthFile' ${config})

echo "scoreratio = ${ratio1}"
echo "pufferfish mapping .. "
echo "/usr/bin/time samtools view -@ 8 -S -b -h <(${pufferfish} align -i ${puffidx} --mate1 ${read1} --mate2 ${read2}  -p 16 -m -o - --scoreRatio ${ratio1}) > ${outdir}/results/pufferfish/${outfile}.bam"
/usr/bin/time samtools view -@ 8 -S -b -h <(${pufferfish} align -i ${puffidx} --mate1 ${read1} --mate2 ${read2}  -p 16 -m -o - --scoreRatio ${ratio1} ) > ${outdir}/results/pufferfish/${outfile}.bam


#/usr/bin/time ${pufferfish} align -i ${puffidx} --mate1 ${read1} --mate2 ${read2}  -p 16 -m -o - --scoreRatio ${ratio1} -k > ${outdir}/results/pufferfish/${outfile}.krkbin


#echo "salmon vb .."
#/usr/bin/time ${salmon} quant -p 16 -la -t ${fasta} -a ${outdir}/results/pufferfish/${outfile}.bam  -o ${outdir}/results/salmon/pufferfish/${outfile}_vb${priorr} --meta --useVBOpt --vbPrior 1e-${priorr}

echo "salmon normal EM .."
echo "/usr/bin/time ${salmon} quant -p 16 -la -t ${fasta} -a ${outdir}/results/pufferfish/${outfile}.bam  -o ${outdir}/results/salmon/pufferfish/${outfile} --meta"
/usr/bin/time ${salmon} quant -p 16 -la -t ${fasta} -a ${outdir}/results/pufferfish/${outfile}.bam  -o ${outdir}/results/salmon/pufferfish/${outfile} --meta

echo "reference level metric results:"
echo "python ${outdir}/scripts/validate_results.py -truer ${truthfile} -est ${outdir}/results/salmon/pufferfish/${outfile}/quant.sf -truthType microbiome"
python ${outdir}/scripts/validate_results.py -truer ${truthfile} -est ${outdir}/results/salmon/pufferfish/${outfile}/quant.sf -truthType microbiome

echo "taxonomy metric results:"
mkdir -p ${outdir}/results/script_results/${outfile}_puff_sal
echo "python ${outdir}/scripts/taxonomicRanks.py -qf ${outdir}/results/salmon/pufferfish/${outfile}/quant.sf -lf ${outdir}/files/rankedlineage.dmp -od ${outdir}/results/script_results/   ${outfile}_puff_sal -dt salmon"
python ${outdir}/scripts/taxonomicRanks.py -qf ${outdir}/results/salmon/pufferfish/${outfile}/quant.sf -lf ${outdir}/files/rankedlineage.dmp -od ${outdir}/results/script_results/${outfile}_puff_sal -dt salmon

echo "python ${outdir}/scripts/validate_results.py -truer ${taxatruthfile} -est ${outdir}/results/script_results/${outfile}_puff_sal/out.tab -truthType taxonomy"
python ${outdir}/scripts/validate_results.py -truer ${taxatruthfile} -est ${outdir}/results/script_results/${outfile}_puff_sal/out.tab -truthType taxonomy




end=`date +%s`
runtime=$((end-start))
printf "\n\nTotal time: $runtime sec\n\n"
