# /PostBusco_treatment (total size = 3.6G)

[1] The Busco raw output are inputed in a homemade python script to perform filtering, quality check, summarize and agregate the data of interest (the orthologous unicopy complete fasta). It's provide 9211 fasta output.
This step also provide the position of each exons in each anotated genes in a metaeuk format. The metaeuk format is converted in gff format using the metaeuktogff.py script and manually sorted by genome position. 
Script exons position in gff format and fasta files are provided.

[2] We reduce the genes dataset to remove genes not present in majority of species and species with a too low Busco score. Leads to a list of 8060 genes and 183 species.
Script to identify genes and species to filter out and list are provided.

[3] We Reduce the 183 species list to keep only one specie per genus (use coverage and availability of their life history traits as a criterion). End up with a list of 144 species (provided). We filter the genes to keep only the sequences from those species.

"X" refers to the gene id

- ### **/script**
 
	- **/Postraitement_busco_V4.py** (input: <busco_split_id>, <species_list>, <path_to_busco_output>, <path_for_outfasta> ; output: <gene_list>, <gene.fasta>, <exon_position>, <busco_summary> )
   	 Postbusco treatment script (step 1). Verify single copy genes, n*3 size of the CDS, count nb of genes per sp and viceversa... 

        - **/gene2sp_filtering.R** (input: <gene2sp_summary>; output: <genelist_filtred_8060genes>, <splist_filtred_183sp>) 
   	R script to identify the genes present in at least 70% of the species and the species containing at least 80% of the remaining genes (step 2).

	- **/metaeuktogff.py** (input: <exon_position>; output: <exon_position.gff>): converte the exons positions in metaeuk format to gff format.

- ### **/data**
	
	- **/9211_no_filter_genes.tar.gz/9211_no_filter_genes/X.fasta** : 9211 genes fasta files, output of the step 1.

	- **/exons_position/exons_position/Xphase1.gff_sorted** : exon position for the 9211 genes, output of the step 1.

	- **/gene2sp_summary** : genexsp presence/absence table for the step 2. Obtained from the merge of the per split gene_list and busco_summary. 

	- **/genelist_filtred_8060genes** : list of 8060 genes contained by at least 70% of species, output of the step 2. 

	- **/splist_filtred_183sp** : list of 183 species presented in at least 80% of the genes, output of the step 2. 

	- **/1genus_144splist** : the 183 species list reduced to 144 species to contain only one specie per gender, output of step 3. 
