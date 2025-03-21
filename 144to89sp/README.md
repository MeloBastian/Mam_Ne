# **/144to89species** (total size = 1.7G)

In an experiment, the 144 species tree was reduced to remove close species using the bayescode cuttree program : /bayescode/data/cuttree <chronogram> <cutoff> <out>.
We choose 0.1 as a cutoff.
When 2 or more species are flagged as too close, one of them is choosen based on its coverage and life history traits availability.
The mapping in then realised as for the 144 species set (see /Mapping).
As there is less species, the genes subset for the mapping run can be larger so the genes are regrouped in 4 * 1500 genes listes rather than 6 * 1000 genes.

- ### **/data**	

	- */90splist_10ma_*:  89 species list
	
	- */1500genes_subset4.conc.treefile_rooted*: rooted treefile from a Iqtree2 run using 1500 genes.
	
	- \*/1500genes_subset*.ali\*: 4 multialignement of 1500 genes. 
	
	- \*/1500genes_subset*.conc\*: 4 concatenation of 1500 genes.
		
	- \*/merge_mam_multigeneglobom.genedsomsuffstat\*: Suffastat from the 4 * 1500 genes mapping, merged in a unique file for fastcoevol.
