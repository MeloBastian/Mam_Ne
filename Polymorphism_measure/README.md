# Polymorphism_measure (total size = 7.4M)

[1] We estimate the number of callable coding sites in the 6002 genes. The inputs of the script comes from different other part of the project.

[2] Compute pS and pN/pS

[3] Study the variability of the pS and pN/pS measure per species by a bootstrap analysis.

- [3.1] Randomly sample with repleacement 6002 genes from the gene list.
   	Create 500 replicates and, for each replicate,  compute a pS and pN/pS.

- [3.2] Aggregation of the results for the different subsets and compute quantiles for pS and pNpS, which are then written in a summary table.

- [3.3] log(pS)~log(pN/pS) graph with bootstrap bars (colors represent the N50 of the genome assemblies).


	
- ### **/script**

	- **/count_coding_callable_6002gene.py** (input: {species_list}, {list_genepersp}, {gene_list}, {bin_gene_mask}; output: {table_nb_coding_callable})
  	Count the number of coding callable positions (step 1).

	- **pnps_6002genes_congruentmasking.py** (input: {table_nb_coding_callable}, {final.vcf}, {species_list}, output: {pNpS_table})
   	pS and pN/pS computation script (step 2).

	- **/table_gene_length_nonalign_aftermasking.py** (input: {genelist}, {splist}, {binarymask_pergene}, {vcf}; output: {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene})
   	Set up the data for the bootstrap analysis. Gives the number of synonymous and non-synonymous sites per genes and the callable size (step 3.1).

	- **/snakefile** for the step 3.1 (input: {genelist}, {specieslist}, {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene}; output: {subset_polym_info} )
   	 Include the scripts:
		
	  - **subsetwriting_redundant.py**
         Write a subset of the gene list by randomly sampling with repeat in the gene list.

	  - **pnps_withmask_fromtable_genexsp.py**
          Adapted script to compute pS and pN/pS only for the genes in the subset.

	- **/bootstrap_calcul.py** (input: {subset_polym_info}, {genes_subset}, output: {bootstrap_summary_table})
   	Merge the 500 subset results and compute summary statistics (pS and pN/pS quantiles) (step 3.2).

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
   	A table with informations like the genome coverage (usefull for step 3.3), n50, number of genes etc 
		
