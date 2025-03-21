# Mam_Ne

This directory contain the key data and script linked to the article _Empirical validation of the nearly neutral theory at divergence and population genomic scale using 150 mammals genomes_ (Bastian, Enard, Lartillot - in prep). 

The aim of this study was to compare estimators of effective population size with measures of selection intensity at phylogenetic and population scales.
To do this, around 150 mammalian genomes were annotated in order to identify approximately 6000 orthologous genes after filtering.
These genes were then aligned, and a phylogeny was reconstructed. They were also used to reconstruct the dN/dS.
In parallel, a vcf file was written for each genome and used to identify any heterozygote position belonging to the genomes studied and their type (synonymous or non-synonymous) in order to calculate a pN/pS and a pS.
Coverage files (bam) were also generated and used to filter the data in order to determine the callable fractions of the genes.
Different filtering steps were performed on the alignments and vcf files.
Alignments, pS,  pN/pS, phylogeny and life history traits were combined and given as input to the Fastcoevol tool, which reconstructed a correlation matrix between traits and the evolution of each trait along the phylogeny.
The correlations are aggregated for graphical representation, with a view to publication.

The "*" in the files names indicate the name of the gene or the specie (depend on context).
The path in the script correspond to the path used for the analysis. You have te replace them to match with your own directory.
The input and output information are not always script argument. The directory of these files is often implemented in the scripts. The input and output indication are only there for a better understanding of how the script works.



-- /PostBusco_treatment 3.6G
(1) The Busco raw output are inputed in a homemade python script to perform filtering, quality check, summarize and agregate the data of interest (the orthologous unicopy complete fasta). It's provide 9211 fasta output.
This step also provide the position of each exons in each anotated genes in a metaeuk format. The metaeuk format is converted in gff format using the metaeuktogff.py script and manually sorted by genome position. 
Script exons position in gff format and fasta files are provided.
(2) We reduce the genes dataset to remove genes not present in majority of species and species with a too low Busco score. Leads to a list of 8060 genes and 183 species.
Script to identify genes and species to filter out and list are provided.
(3) We Reduce the 183 species list to keep only one specie per genus (use coverage and availability of their life history traits as a criterion). End up with a list of 144 species (provided). We filter the genes to keep only the sequences from those species.

	- /script 
		-- /Postraitement_busco_V4.py (input: <busco_split_id>, <species_list>, <path_to_busco_output>, <path_for_outfasta> ; output: <gene_list>, <gene.fasta>, <exon_position>, <busco_summary> )
		  Postbusco treatment script (step 1). Verify single copy genes, n*3 size of the CDS, count nb of genes per sp and viceversa... 
                -- /gene2sp_filtering.R (input: <gene2sp_summary>; output: <genelist_filtred_8060genes>, <splist_filtred_183sp>) 
                  R script to identify the genes present in at least 70% of the species and the species containing at least 80% of the remaining genes (step 2).
		-- /metaeuktogff.py (input: <exon_position>; output: <exon_position.gff>): converte the exons positions in metaeuk format to gff format.
	-- /data
		-- /9211_no_filter_genes/*.fasta : 9211 genes fasta files, output of the step 1.
		-- /exons_position/*phase1.gff_sorted : exon position for the 9211 genes, output of the step 1.
		-- /gene2sp_summary : genexsp presence/absence table for the step 2. Obtained from the merge of the per split gene_list and busco_summary. 
		-- /genelist_filtred_8060genes : list of 8060 genes contained by at least 70% of species, output of the step 2. 
		-- /splist_filtred_183sp : list of 183 species presented in at least 80% of the genes, output of the step 2. 
		-- /1genus_144splist : the 183 species list reduced to 144 species to contain only one specie per gender, output of step 3. 

-- /Genes_alignement 9.3G
(1) We aligned the genes using Prank v.170427 with the command : prank -d={input.fna} -o={output_ali.fna} -DNA -codon. The output are provided
(2) The non-callable positions in the alignement are masked. (the prank results were generated before the obtention of the masking file so we decided to directly masking the prank output rather than reruning it).
(3) The alignements are filtred by differents steps. Only the final filtred alignements are provided in a multialignement (/Mapping/data/mam_subset_*_1000genes.ali).
At the end, there is 6001 genes left.
Different steps:
    (3.1)-- Alignements filtering with HmmCleaner (version 0.180750) but imply first to convert the sequences in amino acide, running Hmmcleaner and then coming back to nucleic sequence using Macse:
        	-- seaview -convert -translate -o {ali_aa.fna} {input_nuc.fna}
        	-- Hmmcleaner {ali_aa.fna}
        	-- macse_v2.06.jar -prog reportMaskAA2NT -align_AA {ali_aa.fna} -align {input_nuc.fna} -mask_AA -
    (3.2)-- Alignements filtering with BMGE (version 1.12): /BMGE.jar -i {input_nuc2.fna} -t CODON -o {output_name}
    (3.3)-- Alignements filtering with a homemade script to remove too short or too gappy files.
    (3.4)-- Alignement analysis with PhylteR version 0.9.7, parameters are given in the logfile (/data/phylter_output)
    (3.5)-- Alignement analysis using Bayescode Version: 0d767fc (https://github.com/bayesiancook/bayescode.git).
	To do so : The data are mapped using the multigeneglobom tools (as described in the /Mapping section). 
      		   The mapping is an input for the tools genebranchdnds (-d {mappingsuffstats}, -t {tree})  which provide a dn/ds measure for each branch and gene in a tab.
      		   This table is analysed to detect outliers branch dn/ds in each genes.
    (3.6)-- Filtering every sequences or genes flagued by Phylter or Bayescode.


	- /script
		-- /mask.py (input : {alignment.fna}, {mask_binary.fna}, {non_alignment.fna}, output : {outname})
		   Masking script (step 2). Masking files cames from /Genome_callable_fraction part.
		-- /ali.py : a package used in mask.py
		-- /alignement_filter_forsnake.py (input: {path_to_genefile}, {gene_name}, output: {out_path})
		   Filter out the too short or too gappy alignements (step 3.3).
  		-- /bayescode_synz_leaf.py (input: {bayescode_output_table}, output: {outliers_sequence_table})
  		   Script to analyse the dn/ds per branch and gene from bayescode (step 3.5). 
		-- /filter_after_phylterandbayescode.py (input: {bayescode_flaged_sequences}, {phylter_flaged_sequences}, {145_splist}, {gene.fasta}; output: {filtred_gene.fasta})
		   Script to filter out the genes or sequences flaged by Phylter or the Bayescode analysis (step 3.6) 
	-- /data
		-- /8060_aligned_fasta/*_filtred_ali.fna : genes after alignement (output of step 1)
                -- /phylter_output : Output of the PhylteR analysis with parameters and flaged sequences.

-- /Mapping 2.8G
After all the alignements filtering step, we endup with a set of 6001 genes.
(1) We agregate them in 6x1000 multialignement file (.ali files) and 6x1000 concatened files (.conc files). This files are provided.
(2) One of the 1000 concatened genes file is used to compute a phylogeny using Iqtree2 and a GTR+4 model (iqtree2 -s {multigene.conc}  -m GTR+G4 -nt 8 -alrt 1000 -B 1000 --boot-trees).
(3) The multialignement and the tree are used for the mapping.
    The mapping is realised using Bayescode Version: 0d767fc (https://github.com/bayesiancook/bayescode.git)
    6 independent mapping were realised corresponding to the 6 set of 1000 genes. The output are then merged in a unique file.

	-- /script 
	There is one script per gene list, they are all similar, only the script for the first list is provided
	        -- /ali_and_conc_afterphylterandbayescode.py (input: {144_splist},{1000_genelist}; output: {multialignement_name.ali}, {multialignement_name.conc})
	           Write a multialignement and a concatenat of 1000 genes. (step 1)
	        -- /mam_multigeneglobom_liste_1 (input: {tree}, {alignement}; output: {out_name})
	        Script for multigenemapping, addapted to be launched on a cluster (step 2)
		-- /mamreadmultigeneglobom_liste_1.slurm (input: {out_mapping_name})
		Treatment of the mappings (step 2).
	-- /data
		-- /list6002genes: final 6001 gene list
		-- /mam_subset_*_1000genes.ali: 6 multialignement of 1000 genes (output of step 1)
		-- /mam_subset_*_1000genes.conc: 6 concatenation of 1000 genes (output of step 1)
		-- 1007forphylo.conc.treefile_rooted: Phylogeny (output of step 2)
		-- /mam_merge_genedsomsuffstat: Suffstat from the 6 * 1000 genes mapping, merged in a unique file for fastcoevol (output of step 3)


-- /144to89species 1.7G
In an experiment, the 144 species tree was reduced to remove close species using the bayescode cuttree program : /bayescode/data/cuttree <chronogram> <cutoff> <out>.
We choose 0.1 as a cutoff.
When 2 or more species are flagged as too close, one of them is choosen based on its coverage and life history traits availability.
The mapping in then realised as for the 144 species set (see /Mapping).
As there is less species, the genes subset for the mapping run can be larger so the genes are regrouped in 4 * 1500 genes listes rather than 6 * 1000 genes.
	-- /data	
		-- /90splist_10ma:  89 species list
		-- /1500genes_subset4.conc.treefile_rooted: rooted treefile from a Iqtree2 run using 1500 genes.
		-- /1500genes_subset*.ali: 4 multialignement of 1500 genes. 
		-- /1500genes_subset*.conc: 4 concatenation of 1500 genes.
		-- /merge_mam_multigeneglobom.genedsomsuffstat: Suffastat from the 4 * 1500 genes mapping, merged in a unique file for fastcoevol.

-- /Genome_callable_fraction 1.4T
(1) The reads are mapped on the genome using bwa_mem2 v2.2.1 and default parameters.
(2) The .bam files obtained from the mapping were used for the variant calling with bcftools v1.13 (mpileup -f {genome_file} {bam_file}) (see /VCF_treatment).
    We obtain 204 .bam files containing the coverage at each position of the genome.
(3) Due to the size of this kind of files (one line per genome position), we resume their information by computing a coverage mean per 100kb windows.
(4) For each resumed bam file, we compute a mean coverage and we identify in a masking file, the 100kb window with a mean coverage 2 times higher or lower than the genome mean.
(5) We intersect the masked genome with the exon positions (from /PostBusco_treatment) to obtain a gene masking.
    We used bcftools subtract -a {exon_position.gff} -b {masked_genome}
(6) We write a binary file for each genes to identify the masked position.
	-- /script
		-- /reducing_coverage_V2.py (input: {deepthfile.gzip}, output: {deepthfile_reduced.gzip})
		   Bam reduction script (step 3).
                -- /mean_med_coveragepergenome.py (input: {deepthfile_reduced.gzip}; output: {summary_file}, {genome_mask.bed})
                  Compute a mean genome coverage and use it to determine masked windows (step 4).
                -- /write_coding_masked_bed_v2.py (input: {144_splist}, {exon_position}, {exon_mask}, output {binary_genemask})
                 rewrite the genes in a binary format corresponding to callable or not position (step 6)
	-- /data 
		-- /depth_files_1.tar.gz and /depth_files_2.tar.gz: Raw genome coverage files (output of step 2).
		-- /binary_genes_masking/*binmask.fasta : binary file (output of step 6).


--  /VCF_treatment 84G
(0) Raw vcf files are generated using bcftools (see /Genome_callable_fraction)
(1) The VCFv4.2 format is not compatible with vcftools (some quality score are in a format not readable for vcftools).
    We reformate the vcf with a homemade script.
(2) We use vcftools to keep only the biallelic SNP : vcftools --vcf {inputvcffile} --max-allele 2 --min-allele 2 --remove-indels --recode --out {outfilename}
(3) We remove the non callable position using the output of the step 4 in the /Genome_callable_fraction part.
    vcftools --vcf {inputvcfile} --exclude-bed {callablefile} --out {outfilename} --recode
(4) We annotate the vcf to characterize each snp as synonymous, non-synonymous or non-coding and filtered out the non-coding positions.
(5) A gene list per species is generated which is needed for next scripts. (provided in /data)
(6) We keep only the snp which belong to the 6002 studied genes and with good quality score.
(7) We count the number of reads per snps and estimate the frequence of each variants and remove outliers snps.

	-- /script
                -- /rewritevcf.py (input: {raw.vcf}; output: {forvcftoold.vcf})
                   Reformate the vcf file for vcftools usage (step 1)
		-- /annotate.bash (input: {vcf}, {genome}, {exon.gff}, {output_directory}, {dir_for_temporary_files})
		   vcf annotation script (step 4).
                -- /annotate_SNPs.py (input: {vcf}, {genome}, {exon.gff} ; output: {annoted.vcf})
                   Python script used in the annotate.bash (step 4).
                -- /formatage.bash (input: {annoted.vcf} {vcf_path} ; output: {annoted_coding.vcf})
                   Remove the non coding SNP (step 4).
                -- /count_seq.py (input: {6002_genes}; output: {list_genepersp})
                   Create a list of gene per species for the 6002 gene list (step 5).
                -- /filteringcodingvcf_vs_homosnp_6002genes.py (input: {sp_list}, {list_genepersp}, {gene_list}, {coding.vcf}; output:{summary_table}, {coding_filtred.vcf})
                   Filter the non homozygote, QUAL<125, GQ <150 and not in the 6002 genes snp (step 6).	
                -- /vcf_allelefreq_table_V2_6002genes.py (input: {gene_list}, {list_genepersp}, {sp_list}, {forvcftools.vcf}, {coding_filtred.vcf}; output: {table_allelicfq_info})
                   Allelic frequency computation (step 7)
		-- /vcf_filtering_fq02to08_6002genes.py (input: {sp_list}, {table_allelicfq_info}, {coding_filtred.vcf}; output: {final.vcf})
		   Allelic frequency filtering (step 7).
		
	-- /data
                -- /vcf_files_lbbe.tar.bz2: Raw vcf files in format VCFv4.2 (output of step 0)
		-- /vcfannoted/*.vcf: vcf annoted, output of the annotate.bash script (output of step 4) 
                -- /table_idgenebysp: list of gene per species for the 6002 gene list (output of step 5)
		-- /finalvcf/Enard_mam_*_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf : final vcf (output of step 7)
		-- /table_sp2lencallcoding_6002genes: summary of the number of coding callable position per species (output of step 7)
		
		
-- Polymorphism_measure 7.4M
(1) We estimate the number of callable coding sites in the 6002 genes. The input of the script cames from different other part of the project.
(2) Compute pS and pN/pS
(3) Study the variability of the pS and pN/pS measure per species by a bootstrap analysis.
    (3.1) Randomly sampling with repeat 6002 genes in the gene list.
    	  Wrote 500 list and compute a pS and pN/pS for each of them.
    (3.2) aggregation of the differents subset results and compute pS and pNpS quantiles in a summary table.
    (3.3) log(pS)~log(pN/pS) graphe with bootstrap bars and n50 coloration
	
	-- /script
		-- /count_coding_callable_6002gene.py (input: {species_list}, {list_genepersp}, {gene_list}, {bin_gene_mask}; output: {table_nb_coding_callable})
		   Count the number of coding callable positions (step 1).
		-- pnps_6002genes_congruentmasking.py (input: {table_nb_coding_callable}, {final.vcf}, {species_list}, output: {pNpS_table})
		   pS and pN/pS computation script (step 2).
		-- /table_gene_length_nonalign_aftermasking.py (input: {genelist}, {splist}, {binarymask_pergene}, {vcf}; output: {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene})
		   Prepare the data for the bootstrap analysis. give the number of synonymous and non-synonymous sites per genes and the callable size (step 3.1)
		-- /snakefile for the step 3.1 (input: {genelist}, {specieslist}, {table_nbcallpospergene}, {table_nbSpospergene}; {table_nbNSpospergene}; output: {subset_polym_info} )
		  Include the scripts:
		  	-- subsetwriting_redundant.py
		  	-- pnps_withmask_fromtable_genexsp.py
		-- /bootstrap_calcul.py (input: {subset_polym_info}, {genes_subset}, output: {bootstrap_summary_table})
	          Merge the subset results and compute summary statistics as the pS and pNpS quantiles (step 3.2)
	        -- /bootstrap_view.R (input: {genomes_summary}, {5%bootstrapquantilesinfo}}
	          plot log(pS) and log(pN/pS) with 5% quantile bootstrap. Color the point in function of the genome coverage and distinguish the 6 species with the suspiciously too light vcf (step 3.3).

	-- data/
		-- Enard_pS_masked_GQ150QUAL125_fq02to08_6002genes: pS and pN/pS per species
		-- /table_genexsp_NScount_6002genes_fq02to08_masked
		  Number of non-synonymous positions per genes (out of step 3.1)
		-- /table_genexsp_Scount_6002genes_fq02to08_masked
		  Number of synonymous positions per genes (out of step 3.1)
		-- /table_spxgenecallablelength_6002genes_fq02to08_masked
		  Number of callable positions per genes (out of step 3.1)
		-- /bootstrap_q5%_psandpnps_6002genes_masked
		   bootstrap summary informations (out of step 3.2)
		-- 144spsummary_genomique_withn50
		  A table with informations like the genome coerage, n50, number of genes etc (usefull for step 3.3)
		

-- /Fastcoevol_analysis 6.5G
(1) Compute a traits table containing Life history traits and polymorphism measures.
(2) Fastcoevol analysis
    It require the suffstats files (/Mapping), a table containing traits of interest at tips (LHT, pS, pN/pS), a rootfile and a tree (/Mapping).
    After a fastcoevol run, the data are readed by a readfastcoevol tool : /bayescode/data/readfastcoevol -x {burning} {everyXpoints} {end} {fastcoevolout}
(3) Commands to transform the .tre output file file in a pdf file containing representing a tree : /coevol/data/drawtree color -th 0.02 -fs 3 file.tre | latex file.tex | dvipdf file.dvi
    
Four analysis were performed with and without the polymorphism of six species and with the 144 species or a subset of 89 more distant species.
The life history traits table was elaged in function of the analyse.

Fastcoevol is a bayescode tools

	-- /script
		-- /thv_Enard.py (input: {pS_pNpS_table}, {anage_lifehistorytraits_database}; output: {traits_table})
		   Agregate LHT and polym (step 1)
		-- /mam_fc_ch1.slurm (input: {traits_table}, {tree}, {suffstats}, {rootfile}, {outputame})
		   Exemple of a Fastcoevol script for a run on a cluster, different chain of analysis were performed, only the script for the first chain is provided (step 2).
	
	-- /data_144sp_allpolym
                -- /rootfile: A file with prior for the traits at the root.
                -- /tableTHV_144mam_masked_6002genes_nogt_lowpsallowed: The table with the life history traits, pS and pN/pS at the tips.
                -- /fastcoevol_output/* , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)
	-- /data_144sp_restrictedpolym
                -- /rootfile: A file with prior for the traits at the root.
                -- /thvwithoutgenerationtime: The table with the life history traits, pS and pN/pS at the tips.
                -- /fastcoevol_output/* , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)
	-- /data_89sp_allpolym
                -- /rootfile: A file with prior for the traits at the root.
                -- /thvwithoutgenerationtime: The table with the life history traits, pS and pN/pS at the tips.
                -- /fastcoevol_output/* , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)
	-- /data-89sp_restrictedpolym	
		-- /rootfile: A file with prior for the traits at the root.
		-- /tableTHV_90mam_masked_6002genes_nogt: The table with the life history traits, pS and pN/pS at the tips.
		-- /fastcoevol_output/* , contains correlation coefficient with posterior probability (.cov), reconstruction of ds and dnds for each branch or nodes and along the tree (.tab and .tre)

-- /Correlation_analysis 65K
(1) Compute dS and dN/dS per branch and leaf for the 89 and 144 species tree. Use the bayescode/get_empirical_branchdnds tools.
(2) Merge the leaf dS and dN/dS with the traits table (/Fastcoevol_analysis).
(2) Independant contrast analysis. 
    Manually report the correlation coefficient and p-value directly in the /ggcorplot_changed.R script. (Report also the fastcoevol results)
(3) Plot the correlation matrice graphe using a modified version of the ggcorrplot R function.
(4) The different correlation matrix figure are manually merged in a unique figure.

	--/script
		-- /pgls.r (input: {traits_dsomenriched_table}, {tree})
		  This script is not automatized. Manual work is needed to choose the 2 traits to correlate and reporting the correlation coefficient and posterior probability.
		  The most important lines are : 
		  	library("caper")
			library(nlme)
		  	IC= comparative.data(phylo, thv, sp, na.omit = F , vcv=TRUE )
		  	modelIC=crunch(log(trait1)~log(trait2), data=IC)
                  	summary(modelIC)
                 -- /ggcorrplot_changed.R
                    This script is not automatized.
                    2 type of matrix are reported : M for the correlation coefficient and PP for the statistical support (posterior probability or p-value).
                    4 version of these matrix are available, corresponding to the 8 studies (2 by matrix).
                    The ggcorrplot function have to be manually modified to accomodate if you analyse a fastcoevol or a pgls analysis (have to modify the sig.level)
	--/data
		-- /tableTHV_90sp_6002genes_lowpsallowed_dsom_forR : LHT+polym+dS+dN/dS table for the 89 species studies with all polymorphism. (step 2)
		-- /tableTHV_90sp_6002genes_dsom_forR: LHT+polym+dS+dN/dS table for the 89 species studies with resticted polymorphism. (step 2)
		-- /tableTHV_144sp_6002genes_lowpsallowed_dsom_forR: LHT+polym+dS+dN/dS table for the 144 species studies with all polymorphism. (step 2)
		-- /tableTHV_144sp_6002genes_dsom_forR: LHT+polym+dS+dN/dS table for the 144 species studies with resticted polymorphism. (step 2)


