#!/bin/bash -l

cellranger_mapping() {

	local sample="$1"

	cellranger count --id="${sample}" \
		--transcriptome=/home/moliva/root_sc_paper/genome/tair10_ref \
		--fastqs=/home/moliva/root_sc_paper/published_datasets/wendrich/fastq/"${sample}"/run1,/home/moliva/root_sc_paper/published_datasets/wendrich/fastq/"${sample}"/run2 \
		--sample="${sample}" \
		--nosecondary \
		--localcores=40 \
		--localmem=300

	
	gunzip -c "${sample}"/outs/filtered_feature_bc_matrix/features.tsv.gz \
	> "${sample}"/outs/filtered_feature_bc_matrix/orig_features.tsv
	awk 'BEGIN {{FS=OFS="\t"}} {{print $1, $1, $3}}' "${sample}"/outs/filtered_feature_bc_matrix/orig_features.tsv \
	| gzip > "${sample}"/outs/filtered_feature_bc_matrix/features.tsv.gz
	rm "${sample}"/outs/filtered_feature_bc_matrix/orig_features.tsv
	
	gunzip -c "${sample}"/outs/raw_feature_bc_matrix/features.tsv.gz \
	> "${sample}"/outs/raw_feature_bc_matrix/orig_features.tsv
	awk 'BEGIN {{FS=OFS="\t"}} {{print $1, $1, $3}}' "${sample}"/outs/raw_feature_bc_matrix/orig_features.tsv \
	| gzip > "${sample}"/outs/raw_feature_bc_matrix/features.tsv.gz
	rm "${sample}"/outs/raw_feature_bc_matrix/orig_features.tsv
	
}

samples=("wendrich_replicate1" \
		"wendrich_replicate2" \
		"wendrich_replicate3")
		
for sample in "${samples[@]}"; do
    cellranger_mapping "$sample"
done