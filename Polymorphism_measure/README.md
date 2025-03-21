# Polymorphism_measure (total size = 7.4M)

[1] We estimate the number of callable coding sites in the 6002 genes. The input of the script cames from different other part of the project.

[2] Compute pS and pN/pS

[3] Study the variability of the pS and pN/pS measure per species by a bootstrap analysis.

   - [3.1] Randomly sampling with repeat 6002 genes in the gene list.
   	Wrote 500 list and compute a pS and pN/pS for each of them.

   - [3.2] Aggregation of the differents subset results and compute pS and pNpS quantiles in a summary table.

   - [3.3] log(pS)~log(pN/pS) graphe with bootstrap bars and n50 coloration


	
- ### **/script**

	- **/count_coding_callable_6002gene.py** (input: {species_list}, {list_genepersp}, {gene_list}, {bin_gene_mask}; output: {table_nb_coding_callable})
  	Count the number of coding callable positions (step 1).

	- **pnps_6002genes_congruentmasking.py** (input: {table_nb_coding_callable}, {final.vcf}, {species_list}, output: {pNpS_table})
   	pS and pN/pS computation script (step 2).

	- **/table_gene_length_nonalign_aftermasking.py** (input: {genelist}, {splist}, {binarymask_pergene}, {vcf}; output: {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene})
   	Prepare the data for the bootstrap analysis. give the number of synonymous and non-synonymous sites per genes and the callable size (step 3.1)

	- **/snakefile** for the step 3.1 (input: {genelist}, {specieslist}, {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene}; output: {subset_polym_info} )
   	 Include the scripts:
		
	  - **subsetwriting_redundant.py**

	  - **pnps_withmask_fromtable_genexsp.py**

	- **/bootstrap_calcul.py** (input: {subset_polym_info}, {genes_subset}, output: {bootstrap_summary_table})
   	Merge the subset results and compute summary statistics as the pS and pNpS quantiles (step 3.2)

	- **/bootstrap_view.R** (input: {genomes_summary}, {5%bootstrapquantilesinfo}}
   		plot log(pS) and log(pN/pS) with 5% quantile bootstrap. Color the point in function of the genome coverage and distinguish the 6 species with the suspiciously too light vcf (step 3.3).

- ### **/data**

	- **Enard_pS_masked_GQ150QUAL125_fq02to08_6002genes**: pS and pN/pS per species

	- **/table_genexsp_NScount_6002genes_fq02to08_masked**
   	Number of non-synonymous positions per genes (out of step 3.1)

	- **/table_genexsp_Scount_6002genes_fq02to08_masked**
   	Number of synonymous positions per genes (out of step 3.1)

	- **/table_spxgenecallablelength_6002genes_fq02to08_masked**
   	Number of callable positions per genes (out of step 3.1)

	- **/bootstrap_q5%_psandpnps_6002genes_masked**
   	bootstrap summary informations (out of step 3.2)

	- **144spsummary_genomique_withn50**
   	A table with informations like the genome coerage, n50, number of genes etc (usefull for step 3.3)
		
