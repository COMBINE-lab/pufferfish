root="/mnt/scratch3/SubmitPuffalign"

##### Aligners #####
puff="/home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish"
bt2="${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2"
star="${root}/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR"	
debga="${root}/tools/deBGA/deBGA"
tools=( "${puff}" "${bt2}" "${star}" "${debga}" )
## Salmon as quantifier ##
salmon="${root}/tools/salmon-latest_linux_x86_64/bin/salmon"

##### Indices #####
puffIdx="${root}/indices/pufferfishindices/viralCG.puffindex"
bt2Idx="${root}/indices/bowtie2indices/viralCG.bowtie2index/bw2"
starIdx="${root}/indices/starindices/viralCG.starindex.fat6/"
debgaIdx="${root}/indices/deBGAindices/viralCG.deBGAindex/"
indices=( "${puffIdx}" "${bt2Idx}" "${starIdx}" "${debgaIdx}" )

##### variable Options #####
options=( 'primary' 'maxNumHits20' 'maxNumHits200' )
puffOpts=( '--bestStrata --primaryAlignment' '--maxNumHits 20' '--maxNumHits200' )
bw2Opts=( '-k 1' '-k 20' '-k 200' )
starOpts=( '--outSAMmultNmax 1' '--outSAMmultNmax 20' '--outSAMmultNmax 200' )
debgaOpts=( '-o 1 -x 1' '-o 20 -x 20' '-o 200 -x 200' )

##### Output directories #####
mkdir -p ${root}/alignments/puffaligner
mkdir -p ${root}/alignments/bowtie2
mkdir -p ${root}/alignments/star
mkdir -p ${root}/alignments/deBGA

##### output statistics #####
if [ -f ${root}/quants/sameSpecies.res ]; then
	rm ${root}/quants/sameSpecies.res
fi

start=`date +%s`
for dataType in 'covid19' 'sars' 'bat2008' ; do
	
	echo ""
	echo "------------------ ${dataType} -------------------"

	dataDir="${root}/data/microbiome_simulation/singleStrainSample/coronas"
	r1="${dataDir}/${dataType}_apr9_1.fastq"
	r2="${dataDir}/${dataType}_apr9_2.fastq"

	for idx in "${!options[@]}"; do
		echo "-----------------------"

		pufName="${root}/alignments/puffaligner/${dataType}.pe.concordant.${options[$idx]}"
		pufCmd="/usr/bin/time ${tools[0]} align --noOrphans --noDiscordant --noDovetail -t 16 ${puffOpts[$idx]} -i ${indices[0]} -1 ${r1} -2 ${r2} -o ${pufName}.sam"

	    bt2Name="${root}/alignments/bowtie2/${dataType}.pe.concordant.${options[$idx]}"
		bt2Cmd="/usr/bin/time ${tools[1]} -t --gbar 1 --no-discordant --no-mixed  -p 16 --ignore-quals --maxins 1000 ${bw2Opts[$idx]} -x ${indices[1]} -1 ${r1} -2 ${r2} -S ${bt2Name}.sam"

    	starName="${root}/alignments/star/${dataType}.pe.concordant.${options[$idx]}"
		mkdir -p ${starName}
		starName="${starName}/star"
		starCmd="/usr/bin/time ${tools[2]}  --runThreadN 16 --alignEndsType EndToEnd --alignIntronMax 1 --outFilterMultimapNmax 1000 --outFilterMismatchNmax 1000 --outFilterMismatchNoverReadLmax 0.20 --alignMatesGapMax 1000 --outMultimapperOrder Random ${starOpts[$idx]} --genomeDir ${indices[2]} --readFilesIn ${r1} ${r2} --outFileNamePrefix ${starName}"

	    debgaName="${root}/alignments/deBGA/${dataType}.pe.concordant.${options[$idx]}"
		debgaCmd="/usr/bin/time ${tools[3]} aln -p 16 -n 1000 -u 1000 -f 0 ${debgaOpts[$idx]} ${indices[3]} ${r1} ${r2} ${debgaName}.sam"
		

		names=( "${pufName}" "${bt2Name}" "${starName}" "${debgaName}" )
		cmds=( "${pufCmd}" "${bt2Cmd}" "${starCmd}" "${debgaCmd}" )
		
		for tidx in "${!tools[@]}"; do
			echo "${tools[$tidx]##*/} ${options[$idx]}"
			cmd="${cmds[$tidx]} > ${names[$tidx]}.log 2>&1"
			echo "${cmd}"
			eval "${cmd}"
			echo "${cmd}" >> ${names[$tidx]}.log
			echo "prep stats .."
			name=${names[$tidx]}
			if [ $tidx -eq 2 ]; then
				name="${name}Aligned.out"
			fi
			echo "${dataType} ${tools[$tidx]##*/} ${options[$idx]}" >> ${root}/quants/sameSpecies.res
			samtools view -@ 16 -F 4 ${name}.sam | cut -f 1 | cut -d'/' -f1 | sort | uniq | wc -l >> ${root}/quants/sameSpecies.res
			samtools view -@ 16 -F 4 ${name}.sam | cut -f 3 | sort | uniq -c >> ${root}/quants/sameSpecies.res
			echo ""
		done
	done
done
echo ""
echo "DONE"
end=`date +%s`
runtime=$((end-start))
echo "Process Time: ${runtime}s"
