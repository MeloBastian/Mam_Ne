### /Genome_callable_fraction (total size = 1.4T)

[1] The reads are mapped on the genome using bwa_mem2 v2.2.1 and default parameters.

[2] The .bam files obtained from the mapping were used for the variant calling with bcftools v1.13 (mpileup -f {genome_file} {bam_file}) (see /VCF_treatment).
   We obtain 204 .bam files containing the coverage at each position of the genome.

[3] Due to the size of this kind of files (one line per genome position), we resume their information by computing a coverage mean per 100kb windows.

[4] For each resumed bam file, we compute a mean coverage and we identify in a masking file, the 100kb window with a mean coverage 2 times higher or lower than the genome mean.

[5] We intersect the masked genome with the exon positions (from /PostBusco_treatment) to obtain a gene masking.
   We used bcftools subtract -a {exon_position.gff} -b {masked_genome}

[6] We write a binary file for each genes to identify the masked position.

"X" refers to the gene id

-# **/script**

	- **/reducing_coverage_V2.py** (input: {deepthfile.gzip}, output: {deepthfile_reduced.gzip})
   		Bam reduction script (step 3).

        - **/mean_med_coveragepergenome.py** (input: {deepthfile_reduced.gzip}; output: {summary_file}, {genome_mask.bed})
   		Compute a mean genome coverage and use it to determine masked windows (step 4).

        - **/write_coding_masked_bed_v2.py** (input: {144_splist}, {exon_position}, {exon_mask}, output {binary_genemask})
   		Rewrite the genes in a binary format corresponding to callable or not position (step 6)

-# **/data** 

	- **/depth_files_1.tar.gz** and **/depth_files_2.tar.gz*: Raw genome coverage files (output of step 2).

	- **/binary_genes_masking/Xbinmask.fasta** : binary file (output of step 6).
