# Mam_Ne

This directory contain the key data and script linked to the article _Empirical validation of the nearly neutral theory at divergence and population genomic scale using 150 mammals genomes_ (Bastian, Enard, Lartillot - in prep). 

**Main goal and analysis**

The aim of this study is to compare estimators of effective population size with measures of selection intensity at phylogenetic (life history traits vs dN/dS) and population (pS vs pN/pS) scales.
To do this, we annotated around 150 mammalian genomes and identify approximately 8000 genes orthologous to mammals (6000 after quality filtering).\
We then align the genes and reconstruct a phylogeny. The genes are also used to reconstruct the dN/dS.\
In parallel, we perform a variant calling analysis for each genome which generate VCF and coverage files.\
The vcf files are used to identify the heterozygote positions in the genes studied and their type (synonymous or non-synonymous) in order to compute a pN/pS and a pS per species.\
The coverage files are used to determine the callable fraction of the data.\
Different filtering steps were performed on the alignments and vcf files to propose a very qualitative dataset.\
Alignments, pS,  pN/pS, phylogeny and life history traits were combined and given as input to the Fastcoevol software, which reconstructed a correlation matrix between traits and the evolution of each trait along the phylogeny.
The correlations are aggregated for graphical representation, with a view to publication.

**Notes for the users**

The different directories corresponds to the differents steps of the analysis. Their is a readme with more detailed informations in each of them.\
In this git repository, only the scripts are provided. The data mentionned in the different readme will be available soon.

The path in the scripts correspond to the path used localy for the analysis. You have te replace them to match with your own directory.
The <input> and <output> information are not always arguments of the script. They are juste indications to the users for a better understanding of how the scripts works.

The pipeline and the logic for using scripts and data are summarised in Figure X.  The colour yellow represents the data (entirely yellow) and script (circled in yellow) that are transmitted. The grey boxes represent the intermediate data.
We have opted not to transmit data that is simple and/or relatively inexpensive to regenerate. All the study scripts are transmitted. If nothing is provided to pass from one file to another, this means that the processing was done by hand, without scripting, in a few simple lines of code.
The coloured backgrounds identify the 8 main stages of the process and refer to the stages presented in the various folders in the directory (except for /144to89sp and /Correlation_analysis).

![](pipeline_scheme.png)


