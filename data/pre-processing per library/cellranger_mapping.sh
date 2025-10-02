#!/bin/bash -l

#This script processes pairs of directories and genome references, utilizing them to execute the 'cellranger count'

cellranger_mapping() {

	local directory="$1"
	local reference="$2"
	
	cd /home/moliva/root_sc_paper/data/raw_data_per_library/"$directory"
	
	# Get paths for fastq subdirectories, separated by a comma
		path_string=""
		# Iterate through subdirectories in the directory
		for subdir in "./fastq"/*; do
			if [ -d "$subdir" ]; then
				path_string+="$(realpath "$subdir"),"
    		fi
		done
		# Remove the trailing comma from the path string
		path_string="${path_string%,}"
	
	# Run cellranger count
		cellranger count --id="cellranger_${directory}" \
		--transcriptome=../../../genome/"${reference}" \
		--fastqs="$path_string" \
		--expect-cells=10000 \
		--nosecondary \
		--localcores=30 \
		--localmem=300

	
	gunzip -c "cellranger_${directory}"/outs/filtered_feature_bc_matrix/features.tsv.gz \
	> "cellranger_${directory}"/outs/filtered_feature_bc_matrix/orig_features.tsv
	awk 'BEGIN {{FS=OFS="\t"}} {{print $1, $1, $3}}' "cellranger_${directory}"/outs/filtered_feature_bc_matrix/orig_features.tsv \
	| gzip > "cellranger_${directory}"/outs/filtered_feature_bc_matrix/features.tsv.gz
	rm "cellranger_${directory}"/outs/filtered_feature_bc_matrix/orig_features.tsv
	
	gunzip -c "cellranger_${directory}"/outs/raw_feature_bc_matrix/features.tsv.gz \
	> "cellranger_${directory}"/outs/raw_feature_bc_matrix/orig_features.tsv
	awk 'BEGIN {{FS=OFS="\t"}} {{print $1, $1, $3}}' "cellranger_${directory}"/outs/raw_feature_bc_matrix/orig_features.tsv \
	| gzip > "cellranger_${directory}"/outs/raw_feature_bc_matrix/features.tsv.gz
	rm "cellranger_${directory}"/outs/raw_feature_bc_matrix/orig_features.tsv
	
}


if [ "$(($# % 2))" -ne 0 ]; then
    echo "This script needs pairs of arguments <directory1> <genomeref1> [<directory2> <genomeref2>...]. Provide an even number of arguments"
    exit 1
fi

# Iterate over pairs of directories and genome references
while [ $# -ge 2 ]; do
    directory="$1"
    reference="$2"
    
    cellranger_mapping "$directory" "$reference"
    
    # Shift the arguments to process the next pair
    shift 2
done
