#  /VCF_treatment (total size = 84G)

[0] Raw vcf files are generated using bcftools (see /Genome_callable_fraction)

[1] The VCFv4.2 format of the vcf is not compatible with vcftools (some quality score are in a format not readable for vcftools).
   We reformate the vcf with a homemade script.

[2] We use vcftools to keep only the biallelic SNP : 
    vcftools --vcf {inputvcffile} --max-allele 2 --min-allele 2 --remove-indels --recode --out {outfilename}

[3] We remove the non callable position using the output of the step 4 in the /Genome_callable_fraction part.
    vcftools --vcf {inputvcfile} --exclude-bed {callablefile} --out {outfilename} --recode

[4] We annotate the vcf to characterize each snp as synonymous, non-synonymous or non-coding and filtered out the non-coding positions.

[5] A gene list per species is generated which is needed for next scripts. 

[6] We keep only the snp which belong to the 6002 studied genes and with good quality score.

[7] We count the number of reads per snps and estimate the frequence of each variants and remove outliers snps.

"X" correspond to the species id

- ### **/script**

  - **/rewritevcf.py** (input: {raw.vcf}; output: {forvcftoold.vcf})
   Reformate the vcf file for vcftools usage (step 1)
	
   - **/annotate.bash** (input: {vcf}, {genome}, {exon.gff}, {output_directory}, {dir_for_temporary_files})
   vcf annotation script (step 4).

   - **/annotate_SNPs.py** (input: {vcf}, {genome}, {exon.gff} ; output: {annoted.vcf})
   Python script used in the annotate.bash (step 4).

    - **/formatage.bash** (input: {annoted.vcf} {vcf_path} ; output: {annoted_coding.vcf})
   Remove the non coding SNP (step 4).

   - **/count_seq.py** (input: {6002_genes}; output: {list_genepersp})
   Create a list of gene per species for the 6002 gene list (step 5).

   - **/filteringcodingvcf_vs_homosnp_6002genes.py** (input: {sp_list}, {list_genepersp}, {gene_list}, {coding.vcf}; output:{summary_table}, {coding_filtred.vcf})
   Filter the non homozygote, QUAL<125, GQ <150 and not in the 6002 genes snp (step 6).	

    - **/vcf_allelefreq_table_V2_6002genes.py** (input: {gene_list}, {list_genepersp}, {sp_list}, {forvcftools.vcf}, {coding_filtred.vcf}; output: {table_allelicfq_info})
   Allelic frequency computation (step 7)

    - **/vcf_filtering_fq02to08_6002genes.py** (input: {sp_list}, {table_allelicfq_info}, {coding_filtred.vcf}; output: {final.vcf})
   Allelic frequency filtering (step 7).
		
- ### **/data**

  - **/vcfannoted.tar.gz/vcfannoted/X.vcf**: vcf annoted, output of the annotate.bash script (output of step 4) 

  - **/table_idgenebysp**: list of gene per species for the 6002 gene list (output of step 5)

  - **finalvcf.tar.gz/finalvcf/Enard_mam_X_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf** : final vcf (output of step 7)

  - **/table_sp2lencallcoding_6002genes** : summary of the number of coding callable position per species (output of step 7)
		
