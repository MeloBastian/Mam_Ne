import sys
import os
import ali

filesp=open("/beegfs/data/mbastian/Enard_postbusco/1genus_144splist_E")
listsp_genus = [line.rstrip("\n") for line in filesp]
aliset = ali.MultiGeneSequenceAlignment(list_name="/beegfs/data/mbastian/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/"+sys.argv[1], dir_name="/beegfs/data/mbastian/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/", format="fasta")
aliset.redefine_taxon_set(listsp_genus, all_taxa=False)
print("before filtering")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

aliset.prune_min_ntaxa(116)

print()
print("after filtering")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

print()
taxon2nbrgene=dict()
for taxon in aliset.taxon_set:
        taxon2nbrgene[taxon]=0
for gene, value in aliset.alignments.items():
        for key in taxon2nbrgene.keys():
                if key in value.get_taxon_set():
                        taxon2nbrgene[key]+=1
print(taxon2nbrgene)

print()
aliset.homogeneize_taxon_sets()

print()
print("after homogeneization")
for gene in aliset.get_gene_list():
    print(gene, aliset.get_gene_ali(gene).get_ntaxa(), aliset.get_gene_ali(gene).get_nsite())

print()

print(len(aliset.alignments))

print("write all alignments to file")
aliset.write_all_to_file(sys.argv[2])
print("concatenate")
concat = aliset.concatenate()
print("total number of sites: ", concat.get_nsite())
print("write concat to file")
concat.write_phylip_to_file(sys.argv[3])



