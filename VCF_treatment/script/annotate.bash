#!/usr/bin/env bash

VCF=$1
GENOME=$2
ANNOT=$3
OUTPUT=$4
#OUTPUT_CODING=$5
DIR_PATH=$5


grep -v '#' $VCF | awk '{print $1,$2,$2+1}' | sed 's/\ /\t/g' > ${DIR_PATH}SNP_pos.bed
#ex : grep -v '#' /beegfs/data/mbastian/VCF_Enard/vcf_D.Enard/Colobus_angolensis.raw.vcf_nonindel.recode.vcf| awk '{print $1,$2,$2+1}' | sed 's/\ /\t/g' >/beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/SNP_pos.bed

cut -f1,3,4,5,6,7 ${ANNOT} > ${DIR_PATH}temp
# cut -f1,3,4,5,6,7 /beegfs/data/mbastian/Enard_postbusco/Refexons/Colobus_angolensis_E_phase1.gff_sorted > /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/temp

/beegfs/home/mbastian/bioinfo/bin/bedtools intersect -a ${DIR_PATH}temp -b ${DIR_PATH}SNP_pos.bed -wa | sed 's/$/\tCDS/' | awk '{print $1,$7,$2,$3,$4,$5,$6}' | sed 's/\ /\t/g' | sort -k3n --unique > ${DIR_PATH}annotation.gff #rm scafold and genes without snp
#/beegfs/home/mbastian/bioinfo/bin/bedtools intersect -a /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/temp -b /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/SNP_pos.bed -wa | sed 's/$/\tCDS/' | awk '{print $1,$7,$2,$3,$4,$5,$6}' | sed 's/\ /\t/g' | sort -k3n --unique > /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/annotation.gff

echo annotation vcf file
/beegfs/home/mbastian/bioinfo/bin/python3 /beegfs/home/mbastian/Scripts/Zoonomia_VCF/annotate_SNPs.py -g $GENOME -vcf $VCF -gff ${DIR_PATH}annotation.gff -o $OUTPUT
#/beegfs/home/mbastian/bioinfo/bin/python3 /beegfs/home/mbastian/Scripts/Zoonomia_VCF/annotate_SNPs.py -g /beegfs/banque/mbastian/Enard/genomes/mam_Colobus_angolensis_E.fna -vcf /beegfs/data/mbastian/VCF_Enard/vcf_D.Enard/Colobus_angolensis.raw.vcf_nonindel.recode.vcf -gff /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/annotation.gff -o /beegfs/data/mbastian/VCF_Enard/vcf_annoted/Colobus_angolensis/Enard_mam_Colobus_angolensis.vcf

#echo remove non coding SNP
## grep NC
#grep -v NC $OUTPUT > ${OUTPUT_CODING}_tmp
#cut -f1,3,4,5,6,7 ${DIR_PATH}annotation.gff > ${DIR_PATH}annotation.gff1 #rm the colone 2 with "\tCDS"
##ajouter colonne position+1 (awk)  awk '{print $1,$2,$2+1,$4,$5,$6,$11,$14,$15}' -> bed format

#echo write bed format
#awk '{print $1, $2, $2+1,$4,$5,$6,$11,$14,$15}' ${OUTPUT_CODING}_tmp | sed 's/\ /\t/g'>${OUTPUT}_coding.bed
## bedtools intersect vcfmodif --> info gene
##rm colonne sans interet gff intriduites dans vcf et colone position +1
##sed ":" info quali
#
#echo keep only interest columns
#/beegfs/home/mbastian/bioinfo/bin/bedtools intersect -b ${OUTPUT}_coding.bed -a ${DIR_PATH}annotation.gff1 -wb | awk '{print $1, $2, $4, $5, $6, $10, $11, $12, $13, $14, $15}'| sed 's/\ /\t/g'  > ${OUTPUT_CODING}
#
#echo tabulation
##tabulate the score
#python3 /beegfs/home/mbastian/Scripts/Zoonomia_VCF/vcf_score_tabulation.py ${OUTPUT_CODING} ${OUTPUT_CODING}_tab
#
#echo all done
#rm annotation.gff
#rm annotation.gff1
#rm temp
#rm ${OUTPUT}_coding.vcf_tmp
#rm SNP_pos.bed

