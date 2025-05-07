#!/bin/bash

rna_path="~/genes_expected_count_DESeq2_H3K27acFormatted_genes_rownames.tsv"
discrete_rna_path="~/genes_expected_count_DESeq2_H3K27acFormatted_dge_binary.tsv"
bigwig_path="~/k27_bigwig_partition0"
segmentation_path="~/segmentation_partition0"
regression_path="~/regression_partition0"
genes_path="~/AllGenes_nonPseudo_nonMit.txt"
gtf_path="~/gencode.v38.annotation.gtf"
chr_sizes="~/hg38_chrSize.txt"

#run expression discretization
Rscript '~/gene_expression_discretization.R' -f "$rna_path"  -g "$gtf_path" -p -b -o "$discrete_rna_path"


mapfile -t gene_ids < "${genes_path}"
echo "Number of filtered genes: ${#gene_ids[@]}"

n=${#gene_ids[@]}
declare -a batch
for (( i=0; i<10000; i++ ))
do
	batch+=(${gene_ids[$i]})
done

N=40
#run STITCHIT on target genes
for gene_id in "${batch[@]}"; do
  ((i=i%N)); ((i++==0)) && wait	
  echo "$gene_id"
  if [ ! -f "${segmentation_path}/Segmentation_${gene_id}_Pearson_10.txt" ]
  then
  	"~/build/core/STITCH" -v 'FALSE' -b "$bigwig_path" -a "$gtf_path" -d "$discrete_rna_path" -o "$rna_path" -s "$chr_sizes" -w 500000 -c 40 -p 0.05 -g "$gene_id" -f "$segmentation_path" -z 10 -r 1100000 -t 5000 -m "Pearson" & 
  fi
done

Rscript "~/Two_level_learning_elasticNet.R" --dataDir="$segmentation_path" --outDir="$regression_path" --response="Expression" --cores=40 --alpha=0.05 --testsize=0.2 --regularisation='E' --innerCV=6 --outerCV=6 --performance='TRUE' --leaveOneOutCV='FALSE' --asRData='TRUE' --randomise='FALSE' --logResponse='TRUE' --ftest='FALSE' --coefP=1
