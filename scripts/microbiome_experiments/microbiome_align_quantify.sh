root="/mnt/scratch3/SubmitPuffalign"

mkdir -p ${root}/alignments/puffaligner
mkdir -p ${root}/alignments/bowtie2
mkdir -p ${root}/alignments/star

#names=( 'primary' 'primary.allowDove' 'maxNumHits20' 'maxNumHits20.allowDove' )
#puffOpts=( '--noDovetail --bestStrata --primaryAlignment' '--bestStrata --primaryAlignment' '--noDovetail --maxNumHits 20' '--maxNumHits 20' )
#bw2Opts=( '-k 1' '-k 1 --dovetail' '-k 20' '-k 20 --dovetail' )


names=( 'primary' 'maxNumHits20')
puffOpts=( '--bestStrata --primaryAlignment' '--maxNumHits 20' )
bw2Opts=( '-k 1' '-k 20' )
starOpts=( '--outSAMmultNmax 1' '--outSAMmultNmax 20' )
#names=( 'maxNumHits20.genomicReads' )
#puffOpts=( '--maxNumHits 20 --genomicReads' )

echo "Distinguishing Similar Sequences"
echo ""
echo ""
singleStrainDir="${root}/data/microbiome_simulation/singleStrainSample/coronas"
for dataType in 'covid19' 'sars' 'bat2008' ; do
	echo "------------------ ${dataType} -------------------"
	echo ""

	for idx in "${!names[@]}"; do
		echo "${names[$idx]}"
		echo ""

		echo "Puffaligner"
		name="${root}/alignments/puffaligner/${dataType}.pe.concordant.${names[$idx]}"
		cmd="/usr/bin/time /home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish align -i ${root}/indices/pufferfishindices/viralCG.puffindex/ -1 ${singleStrainDir}/${dataType}_apr9_1.fastq -2 ${singleStrainDir}/${dataType}_apr9_2.fastq -o ${name}.sam  --noOrphans --noDiscordant --noDovetail ${puffOpts[$idx]} -t 16 > ${name}.log 2>&1"
		echo ${cmd}
		eval ${cmd}
		echo ${cmd} >> ${name}.log
		echo ""
	

		echo "Bowtie2"
		name="${root}/alignments/bowtie2/${dataType}.pe.concordant.${names[$idx]}"
		cmd="/usr/bin/time ${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2 ${bw2Opts[$idx]} --gbar 1 --no-discordant --no-mixed  -x ${root}/indices/bowtie2indices/viralCG.bowtie2index/bw2 -1 ${singleStrainDir}/${dataType}_apr9_1.fastq -2 ${singleStrainDir}/${dataType}_apr9_2.fastq -S ${name}.sam -p 16 --ignore-quals --maxins 1000 > ${name}.log 2>&1"
		echo ${cmd}
		eval ${cmd}
		echo ${cmd} >> ${name}.log
		echo ""

		echo "STAR"
		name="${root}/alignments/star/${dataType}.pe.concordant.${names[$idx]}"
		mkdir -p ${name}
		cmd="/usr/bin/time ${root}/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir ${root}/indices/starindices/viralCG.starindex.fat6/ --readFilesIn ${singleStrainDir}/${dataType}_apr9_1.fastq ${singleStrainDir}/${dataType}_apr9_2.fastq --runThreadN 16 --alignEndsType EndToEnd --outFilterMultimapNmax 1000 --outFilterMismatchNmax 1000 --outFilterMismatchNoverReadLmax 0.12 --alignMatesGapMax 1000 ${starOpts[$idx]} --outFileNamePrefix ${name}/  > ${name}/log 2>&1"
		echo ${cmd}
		eval ${cmd}
		echo ${cmd} >> ${name}/log
		echo ""

		echo "-------------------"
		echo ""

	done

done

############################################################################################################################################################################
echo "Complex Metagenomic Sample Alignment"
echo ""
echo ""


dataType="soil_SRR1094822"
dataType="jiandongPeninsula_SRR11283975"
dataType="uzonCaldera_SRR11496426"

dataDir="${root}/data/microbiome_simulation/bulkSample/${dataType}"
r1="${dataDir}/simulatedReads_1.fastq"
r2="${dataDir}/simulatedReads_2.fastq"

for idx in "${!names[@]}"; do
	echo "${names[$idx]}"
	echo ""

	echo "Puffaligner"
	name="${root}/alignments/puffaligner/${dataType}.pe.concordant.${names[$idx]}"
	cmd="/usr/bin/time /home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish align -i ${root}/indices/pufferfishindices/microbiome_k31.resubmission.puffindex.sparse/ -1 ${r1} -2 ${r2} -o ${name}.sam  --noOrphans --noDiscordant --noDovetail ${puffOpts[$idx]} -t 16 > ${name}.log 2>&1"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${name}.log
	samtools view -@ 16 -bS ${name}.sam > ${name}.bam
	rm ${name}.sam
	echo ""
########## salmon ##########
	echo "Salmon"
	quantname="${root}/quants/salmon/${dataType}.puff.${names[$idx]}"
	mkdir -p ${quantname}
	cmd="/usr/bin/time ${root}/tools/salmon-latest_linux_x86_64/bin/salmon quant -p 16 -la -a ${name}.bam -t ${root}/microbiome/bacteria_virus/library/bacteria/selected_complete_genome.fna -o ${quantname} 2>${quantname}/log"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${quantname}/log



	echo ""
	echo "Bowtie2"
	name="${root}/alignments/bowtie2/${dataType}.pe.concordant.${names[$idx]}"
	cmd="/usr/bin/time ${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2 ${bw2Opts[$idx]} --gbar 1 --no-discordant --no-mixed  -x ${root}/indices/bowtie2indices/microbiome.bowtie2index/bw2 -1 ${r1} -2 ${r2} -S ${name}.sam -p 16 --ignore-quals --maxins 1000 > ${name}.log 2>&1"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${name}.log
	echo ""

######### salmon ##########
	echo "Salmon"
	cat <(samtools view -H ${name}.sam) <(samtools view ${name}.sam | awk -v OFS='\t' '{gsub(/\/1|\/2/,"",$1); print $0;}') | samtools view -@ 16 -bS > ${name}.bam
	rm ${name}.sam
	quantname="${root}/quants/salmon/${dataType}.bt2.${names[$idx]}"
	mkdir -p ${quantname}
	cmd="/usr/bin/time ${root}/tools/salmon-latest_linux_x86_64/bin/salmon quant -p 16 -la -a ${name}.bam -t ${root}/microbiome/bacteria_virus/library/bacteria/selected_complete_genome.fna -o ${quantname} 2>${quantname}/log"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${quantname}/log


	echo "-------------------"
	echo ""

done

