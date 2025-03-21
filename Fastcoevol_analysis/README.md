
# **/Fastcoevol_analysis** (total size = 6.5G)

[1] Compute a traits table containing Life history traits and polymorphism measures.

[2] Fastcoevol analysis
   It require the suffstats files (/Mapping), a table containing traits of interest at tips (LHT, pS, pN/pS), a rootfile and a tree (/Mapping).
   After a fastcoevol run, the data are readed by a readfastcoevol tool : /bayescode/data/readfastcoevol -x {burning} {everyXpoints} {end} {fastcoevolout}

[3] Commands to transform the .tre output file file in a pdf file containing representing a tree : /coevol/data/drawtree color -th 0.02 -fs 3 file.tre | latex file.tex | dvipdf file.dvi

Four analysis were performed with and without the polymorphism of six species and with the 144 species or a subset of 89 more distant species.
The life history traits table was elaged in function of the analyse.

Fastcoevol is a bayescode tools

- ### **/script**

	- **/thv_Enard.py** (input: {pS_pNpS_table}, {anage_lifehistorytraits_database}; output: {traits_table})
   		Agregate LHT and polym (step 1)
		
	-  **/mam_fc_ch1.slurm** (input: {traits_table}, {tree}, {suffstats}, {rootfile}, {outputame})
   		Exemple of a Fastcoevol script for a run on a cluster, different chain of analysis were performed, only the script for the first chain is provided (step 2).
	
- ### **/data_144sp_allpolym**

	- **/rootfile**: A file with prior for the traits at the root.

	- **/tableTHV_144mam_masked_6002genes_nogt_lowpsallowed**: The table with the life history traits, pS and pN/pS at the tips.

	- **/fastcoevol_output/** , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)
	
- ### **/data_144sp_restrictedpolym**

        - **/rootfile**: A file with prior for the traits at the root.

        - **/thvwithoutgenerationtime**: The table with the life history traits, pS and pN/pS at the tips.

        - **/fastcoevol_output/** , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)

- ###  **/data_89sp_allpolym**

        - **/rootfile**: A file with prior for the traits at the root.

        - **/thvwithoutgenerationtime**: The table with the life history traits, pS and pN/pS at the tips.

        - **/fastcoevol_output/** , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)

- ### **/data-89sp_restrictedpolym**	

	- **/rootfile**: A file with prior for the traits at the root.

	- **/tableTHV_90mam_masked_6002genes_nogt**: The table with the life history traits, pS and pN/pS at the tips.

	- **/fastcoevol_output/** , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)
