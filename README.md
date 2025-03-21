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


add the pipeline figure
