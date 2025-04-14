# /Genes_alignement (total size = 9.3G)

[1] We aligned the genes using Prank v.170427 with the command : prank -d={input.fna} -o={output_ali.fna} -DNA -codon. 

[2] The non-callable positions in the alignement are masked. (the prank results were generated before the obtention of the masking file so we decided to directly mask the prank output rather than reruning it).

[3] The alignements are filtered by differents steps.
At the end, there is 6001 genes left.

Different alignment filtering steps:

   	[3.1] Alignements filtering with HmmCleaner (version 0.180750) but imply first to convert the sequences in amino acide, running Hmmcleaner and then coming back to nucleic sequence using Macse:
    		- seaview -convert -translate -o {ali_aa.fna} {input_nuc.fna}
        	- Hmmcleaner {ali_aa.fna}
        	- macse_v2.06.jar -prog reportMaskAA2NT -align_AA {ali_aa.fna} -align {input_nuc.fna} -mask_AA -
 
   	[3.2] Alignements filtering with BMGE (version 1.12): /BMGE.jar -i {input_nuc2.fna} -t CODON -o {output_name}

   	[3.3] Alignements filtering with a homemade script to remove too short or too gappy files (alignement_filter_forsnake.py).

   	[3.4] Alignement analysis with PhylteR version 0.9.7, parameters are given in the logfile (/data/phylter_output)

   	[3.5] Alignement analysis using Bayescode Version: 0d767fc (https://github.com/bayesiancook/bayescode.git).
 	      To do so :
                -  The data are mapped using the multigeneglobom tools (as described in the /Mapping section). 
                -  The mapping is an input for the tools genebranchdnds (-d {mappingsuffstats}, -t {tree})  which provide a dN/dS measure for each branch and gene in a tab.
                -  This table is analysed to detect outliers branch dN/dS in each genes.

   	[3.6] Filtering every sequences or genes flagged by Phylter or Bayescode.

"X" correspond to the genes id

- ### **/script**

	- **/mask.py** (input : {alignment.fna}, {mask_binary.fna}, {non_alignment.fna}, output : {outname})
   	Masking script (step 2). Masking files cames from the /Genome_callable_fraction part.

	- **/ali.py** : a package used in mask.py, developped by Nicolas Lartillot.

	- **/alignement_filter_forsnake.py** (input: {path_to_genefile}, {gene_name}, output: {out_path})
   	Filter out the too short or too gappy alignements (step 3.3).

	- **/bayescode_synz_leaf.py** (input: {bayescode_output_table}, output: {outliers_sequence_table})
  	Script to analyse the dN/dS per branch and gene from bayescode (step 3.5). 

	- **/filter_after_phylterandbayescode.py** (input: {bayescode_flaged_sequences}, {phylter_flaged_sequences}, {145_splist}, {gene.fasta}; output: {filtred_gene.fasta})
  	Script to filter out the genes or sequences flaged by Phylter or the Bayescode analysis (step 3.6) 

- ### **/data**

	- **/8060_aligned_fasta.tar.gz/8060_aligned_fasta/X_filtred_ali.fna** : genes after alignement (output of step 1)

	- **/phylter_output** : Output of the PhylteR analysis with parameters and flaged sequences.
