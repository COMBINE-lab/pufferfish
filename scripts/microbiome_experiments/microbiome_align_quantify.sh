root="/mnt/scratch3/SubmitPuffalign"

mkdir -p ${root}/alignments/puffaligner
mkdir -p ${root}/alignments/bowtie2

#names=( 'primary' 'primary.allowDove' 'maxNumHits20' 'maxNumHits20.allowDove' )
#puffOpts=( '--noDovetail --bestStrata --primaryAlignment' '--bestStrata --primaryAlignment' '--noDovetail --maxNumHits 20' '--maxNumHits 20' )
#bw2Opts=( '-k 1' '-k 1 --dovetail' '-k 20' '-k 20 --dovetail' )


names=( 'primary' 'maxNumHits20')
puffOpts=( '--bestStrata --primaryAlignment' '--maxNumHits 20' )
bw2Opts=( '-k 1' '-k 20' )

echo "Distinguishing Similar Sequences"
echo ""
echo ""
for dataType in 'covid19' 'sars' 'bat2008' ; do
	echo "------------------ ${dataType} -------------------"
	echo ""

	for idx in "${!names[@]}"; do
		echo "${names[$idx]}"
		echo ""

		echo "Puffaligner"
		name="${root}/alignments/puffaligner/${dataType}.pe.concordant.${names[$idx]}"
		cmd="/usr/bin/time /home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish align -i ${root}/indices/pufferfishindices/viralCG.puffindex/ -1 ${root}/data/microbiome_simulation/${dataType}_apr9_1.fastq -2 ${root}/data/microbiome_simulation/${dataType}_apr9_2.fastq -o ${name}.sam  --noOrphans --noDiscordant --noDovetail ${puffOpts[$idx]} -t 16 > ${name}.log 2>&1"
		echo ${cmd}
		eval ${cmd}
		echo ${cmd} >> ${name}.log
		echo ""
	

		echo "Bowtie2"
		name="${root}/alignments/bowtie2/${dataType}.pe.concordant.${names[$idx]}"
		cmd="/usr/bin/time ${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2 ${bw2Opts[$idx]} --gbar 1 --no-discordant --no-mixed  -x ${root}/indices/bowtie2indices/viralCG.bowtie2index/bw2 -1 ${root}/data/microbiome_simulation/${dataType}_apr9_1.fastq -2 ${root}/data/microbiome_simulation/${dataType}_apr9_2.fastq -S ${name}.sam -p 16 --ignore-quals --maxins 1000 > ${name}.log 2>&1"
		echo ${cmd}
		eval ${cmd}
		echo ${cmd} >> ${name}.log

		echo "-------------------"
		echo ""

	done

done

############################################################################################################################################################################
echo "Complex Metagenomic Sample Alignment"
echo ""
echo ""

dataDir="${root}/data/microbiome_simulation/jiandongPeninsula_SRR11283975"
r1="${dataDir}/simulatedReads_1.fastq"
r2="${dataDir}/simulatedReads_2.fastq"
dataType="jiandongPeninsula_SRR11283975"

for idx in "${!names[@]}"; do
	echo "${names[$idx]}"
	echo ""

	echo "Puffaligner"
	name="${root}/alignments/puffaligner/${dataType}.pe.concordant.${names[$idx]}"
	cmd="/usr/bin/time /home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish align -i ${root}/indices/pufferfishindices/microbiome_k31.resubmission.puffindex.sparse/ -1 ${r1} -2 ${r2} -o ${name}.sam  --noOrphans --noDiscordant --noDovetail ${puffOpts[$idx]} -t 16 > ${name}.log 2>&1"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${name}.log
	echo ""
########## salmon ##########
	echo "Salmon"
	quantname="${root}/quants/salmon/${dataType}.puff.${names[$idx]}"
	mkdir -p ${quantname}
	cmd="/usr/bin/time ${root}/tools/salmon-latest_linux_x86_64/bin/salmon quant -p 16 -la -a ${name}.sam -t ${root}/microbiome/bacteria_virus/library/bacteria/selected_complete_genome.fna -o ${quantname} 2>${quantname}/log"
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
	cat <(samtools view -H ${name}.sam) <(samtools view ${name}.sam | awk -v OFS='\t' '{gsub(/\/1|\/2/,"",$1); print $0;}') > ${name}_2.sam
	quantname="${root}/quants/salmon/${dataType}.bt2.${names[$idx]}"
	mkdir -p ${quantname}
	cmd="/usr/bin/time ${root}/tools/salmon-latest_linux_x86_64/bin/salmon quant -p 16 -la -a ${name}_2.sam -t ${root}/microbiome/bacteria_virus/library/bacteria/selected_complete_genome.fna -o ${quantname} 2>${quantname}/log"
	echo ${cmd}
	eval ${cmd}
	echo ${cmd} >> ${quantname}/log


	echo "-------------------"
	echo ""

done


