#!/bin/bash -l

#This script processes pairs of library directories and ecotype SNP information, utilizing them to count SNPS per cell using sctools


run_sctools() {

	local directory="$1"
	local ecotype="$2"
	
	cd /home/moliva/root_sc_paper/data/raw_data_per_library/"$directory"
	
	mkdir genotyping
	
	sctools countsnps -b "cellranger_${directory}"/outs/possorted_genome_bam.bam \
	-s /home/moliva/root_sc_paper/genome/"${ecotype}_snp_last_500bp.tsv" \
	-o genotyping/"${directory}_snps.tsv" \
	-p 20
	
}

if [ "$(($# % 2))" -ne 0 ]; then
    echo "This script needs pairs of arguments. Provide an even number of arguments"
    exit 1
fi

# Iterate over pairs of directories and ecotype SNP files
while [ $# -ge 2 ]; do
    directory="$1"
    ecotype="$2"
    
    run_sctools "$directory" "$ecotype"
    
    # Shift the arguments to process the next pair
    shift 2
done
