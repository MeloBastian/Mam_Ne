#!/usr/bin/env bash


OUTPUT=$1
OUTPUT_CODING=$2
DIR_PATH=$3
#OUTPUT_CODING_tab=$4

echo remove non coding SNP
## grep NC
awk '{if ($6 != "NC") print}' $OUTPUT > ${OUTPUT_CODING}_tmp
cut -f1,3,4,5,6,7 ${DIR_PATH}/annotation.gff > ${DIR_PATH}/annotation.gff1 #rm the colone 2 with "\tCDS"
##ajouter colonne position+1 (awk)  awk '{print $1,$2,$2+1,$4,$5,$6,$11,$14,$15}' -> bed format

echo write bed format
awk '{print $1, $2, $2+1,$4,$5,$6,$11,$14,$15}' ${OUTPUT_CODING}_tmp | sed 's/\ /\t/g'>${OUTPUT}_coding.bed


echo keep only interest columns
/beegfs/home/mbastian/bioinfo/bin/bedtools intersect -b ${OUTPUT}_coding.bed -a ${DIR_PATH}/annotation.gff1 -wb | awk '{print $1, $2, $4, $5, $6, $10, $11, $12, $13, $14, $15}'| sed 's/\ /\t/g'  > ${OUTPUT_CODING}

