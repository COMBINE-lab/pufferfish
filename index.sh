#!/bin/bash
#!/usr/bin/env

echo "Building Pufferfish Index"
echo "Pipeline is fixFasta -> twoPaco -> twoPaCo dump -> pufferize -> pufferfish index"
echo ""

PUFFERFISH=$(jq -r '.pufferfish' config.json)/build/src
TWOPACO=$(jq -r '.twopaco' config.json)/build/graphconstructor/twopaco
TWOPACO_DUMP=$(jq -r '.twopaco' config.json)/build/graphdump/graphdump

echo "pufferfish: $PUFFERFISH"
echo "twopaco: $TWOPACO"

K=$(jq -r '.ksize' config.json)
OUTPUT_DIR=$(jq -r '.output_dir' config.json)
INPUT_FILES=$(jq -r '.input.input_fasta' config.json)
TMP=$(jq -r '.tmp_dir' config.json)

echo "k: $K"
echo "output directory: $OUTPUT_DIR"
echo "input directory/file: $INPUT_FILES"
echo "tmp directory: $TMP"

THREAD=$(jq -r '.num_of_threads' config.json)
echo "num of threads: $THREAD"

IS_INPUT_DIRECTORY=$(jq -r '.input.is_input_a_directory_to_fasta_files' config.json)

mkdir -p $OUTPUT_DIR
mkdir -p $TMP

if [ $IS_INPUT_DIRECTORY = true ]; then
	if [ ! -d $INPUT_FILES ]; then
		echo "ERROR: $INPUT_FILES does not exist"
		exit 1
	fi
elif [ $IS_INPUT_DIRECTORY = false ]; then
    if [ ! -f $INPUT_FILES ]; then
	    echo "ERROR: $INPUT_FILES does not exist"
		exit 1
	fi
fi
if [ ! -d $PUFFERFISH ]; then                                                           
    echo "ERROR: $PUFFERFISH does not exist. Make sure you've compiled pufferfish"                                               
    exit 1                                                                               
fi

if [ ! -f $TWOPACO ]; then                                                        
    echo "ERROR: $TWOPACO does not exist. Make sure you've compiled TwoPaCo"                                               
    exit 1                                 
fi                                                                                      

printf "\nStarting the pipeline\n"

printf "\nFixFasta:\n"
bname="ref"
if [ $IS_INPUT_DIRECTORY = true ]; then
	echo "Input is the path to the directory of fasta files"
	wc -l <(cat `ls -d $INPUT_FILES/*.*`)
	bname=$bname"_fixed"
	$PUFFERFISH/fixFasta --klen $K --input <(cat `ls -d $INPUT_FILES/*.*`) --output $OUTPUT_DIR/${bname}.fa
elif [ $IS_INPUT_DIRECTORY = false ]; then
	echo "Input is the path to the concatenated fasta file"
	bname=`basename ${INPUT_FILES} | sed 's/\.[^\.]*$//'`
	bname=$bname"_fixed"
	echo "$bname"
	echo "$PUFFERFISH/fixFasta --klen $K --input $INPUT_FILES --output $OUTPUT_DIR/${bname}.fa"
	$PUFFERFISH/fixFasta --klen $K --input $INPUT_FILES --output $OUTPUT_DIR/${bname}_fixed.fa
else
	echo "ERROR: Wrong format for a boolean value in json file config.json"
	exit 1
fi


printf "\nTwoPaCo Junction Detection:\n"
echo "$TWOPACO -k $K -t $THREAD -f 20 "$OUTPUT_DIR/$bname.fa" --outfile $OUTPUT_DIR/$bname"_dbg.bin" --tmpdir $TMP"
/usr/bin/time $TWOPACO -k $K -t $THREAD -f 20 "$OUTPUT_DIR/$bname.fa" --outfile $OUTPUT_DIR/$bname"_dbg.bin" --tmpdir $TMP 

printf "\nTwoPaCo Dump:\n"
echo "$TWOPACO_DUMP -k $K -s "$OUTPUT_DIR/$bname.fa" -f gfa1 $OUTPUT_DIR/$bname"_dbg.bin" > $OUTPUT_DIR/$bname.gfa"
/usr/bin/time $TWOPACO_DUMP -k $K -s "$OUTPUT_DIR/$bname.fa" -f gfa1 $OUTPUT_DIR/$bname"_dbg.bin" > $OUTPUT_DIR/$bname.gfa

printf "\nPufferfish Puffferize:\n"
printf "$PUFFERFISH/pufferize -k $K -g $OUTPUT_DIR/$bname.gfa -f $OUTPUT_DIR/$bname.fa -o $OUTPUT_DIR/$bname\_pufferized.gfa" 
/usr/bin/time $PUFFERFISH/pufferize -k $K -g "$OUTPUT_DIR/$bname.gfa" -f "$OUTPUT_DIR/$bname.fa" -o $OUTPUT_DIR/$bname"_pufferized.gfa"

printf "\n\nPufferfish Build Index:\n"
PUFFERFISH_INDEX="$OUTPUT_DIR/$bname.puffidx"
mkdir -p $PUFFERFISH_INDEX
echo "$PUFFERFISH/pufferfish index -k $K -g $OUTPUT_DIR/$bname"_pufferized.gfa" -r "$OUTPUT_DIR/$bname.fa" -o $PUFFERFISH_INDEX"
/usr/bin/time $PUFFERFISH/pufferfish index -k $K -g $OUTPUT_DIR/$bname"_pufferized.gfa" -r "$OUTPUT_DIR/$bname.fa" -o $PUFFERFISH_INDEX
