#! /usr/bin/python3
import sys
import os
import random
import re

univ_trans = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_codon(codon):
    if codon in univ_trans:
        return univ_trans[codon]
    else:
        return "?"

def translate_nucseq(nucseq):
    nsite = len(nucseq)
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    return "".join([translate_codon(nucseq[3*i:3*(i+1)]) for i in range(naa)])

def get_gc3(nucseq):
    nsite = len(nucseq)
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    ngc = 0
    ntot = 0
    for i in range(naa):
        if nucseq[3*i+2] in "gcGC":
            ngc += 1
        if nucseq[3*i+2] in "atgcATGC":
            ntot += 1
    f = 0
    if ntot:
        f = ngc / ntot
    return f

class SequenceAlignment:

    def __init__(self, file_name = "", file_stream = "", ali = "", format = "phylip", fromdict = dict()):
        self.nsite = 0
        self.ali = dict()
        if file_name:
            if format == "phylip":
                self.read_sequential_phylip_from_file(file_name)
            if format == "fasta":
                self.read_fasta_from_file(file_name)
            if format == "nexus":
                self.read_sequential_nexus_from_file(file_name)

        if file_stream:
            if format != "phylip":
                print("error in SequenceAlignment: stream constructor only accepts phylip format")
                sys.exit()
            self.read_sequential_phylip_from_stream(file_stream)

        if len(fromdict):
            self.ali = fromdict
            self.nsite = max([len(seq) for tax,seq in self.ali.items()])

    def get_taxon_set(self):
        return {tax for (tax,seq) in self.ali.items()}

    def get_ntaxa(self):
        return len(self.ali)

    def get_total_size(self):
        return sum([len(seq) for (tax,seq) in self.ali.items()])

    def get_nsite(self):
        return self.nsite

    def get_seq(self, tax):
        if tax not in self.ali:
            print("error in SequenceAlignment.get_seq: taxon not in alignment")
            sys.exit(1)
        return self.ali[tax]

    def read_sequential_phylip_from_file(self, filename):
        with open(filename, 'r') as infile:
            self.read_sequential_phylip_from_stream(infile)

    def read_sequential_phylip_from_stream(self, infile):
        (ntax,npos) = infile.readline().rstrip('\n').split()
        ntaxa = int(ntax)
        self.nsite = int(npos)

        self.ali = dict()
        for i in range(ntaxa):
            line = infile.readline()
            if not line:
                print("error when reading phylip")
                sys.exit()
            pair = line.rstrip('\n').split()
            if self.nsite:
                if len(pair) != 2:
                    print("error in readphylip: should have 2 fields per line")
                    print(pair)
                    sys.exit()
                if len(pair[1]) != self.nsite:
                    print("error in read_phylip: non matching number of sites")
                    sys.exit()
                self.ali[pair[0]] = pair[1]
            else:
                self.ali[pair[0]] = ""

    def read_sequential_nexus_from_file(self, filename):
        with open(filename, 'r') as infile:
            self.read_sequential_nexus_from_stream(infile)

    def read_sequential_nexus_from_stream(self, infile):
        infile.readline()
        infile.readline()
        line = infile.readline().rstrip('\n')
        match_header = r"^\s+dimensions\s+ntax=(\d+)\s+nchar=(\d+);$"
        m = re.match(match_header, line)
        if not m:
            print("error: no match for nexus")
            print(line)
        ntaxa = int(m.group(1))
        self.nsite = int(m.group(2))

        infile.readline()
        infile.readline()

        self.ali = dict()
        for i in range(ntaxa):
            line = infile.readline()
            if not line:
                print("error when reading phylip")
                sys.exit()
            pair = line.rstrip('\n').split()
            if len(pair) != 2:
                print("error in readphylip: should have 2 fields per line")
                print(pair)
                sys.exit()
            if len(pair[1]) != self.nsite:
                print("error in read_phylip: non matching number of sites")
                print(pair[1])
                sys.exit()
            self.ali[pair[0]] = pair[1]

    def mask_ali(self, masks, unali):
        print("taxon\tunali_length\tali_length\t#non_match")
        for (mtax, mask) in masks.ali.items():
            tax = mtax + "_E"
            if tax not in unali.ali:
                print("error: did not find taxon {0} in unaligned sequences".format(tax)) 
                sys.exit()
            unaliseq = unali.ali[tax]
            unaliseq = unaliseq.upper()
            if tax not in self.ali:
                print("error: did not find taxon {0} in alignement".format(tax))
                sys.exit()
            seq = self.ali[tax]
            seq_n = len(seq)
            new_seq = ""
            mask_n = len(mask)

            if mask_n != len(unaliseq):
                print("mask and unaligned sequence have different lengths: {0}/{1}".format(mask_n, len(unaliseq)))
                sys.exit()

            seq_i = 0
            mask_i = 0
            nonmatch = 0;
            while (seq_i < seq_n):
                if seq[seq_i] == '-':
                    new_seq += '-'
                else:
                    if mask_i >= mask_n:
                        print("error: mask length overflow")
                        sys.exit()

                    if unaliseq[mask_i] != seq[seq_i]:
                        nonmatch += 1
                        if seq[seq_i] != 'N':
                            print("error: non matching characters between position {0} of unaligned sequence (char {1}) and position {2} of aligned sequence (char {3})".format(mask_i, unaliseq[mask_i], seq_i, seq[seq_i]))
                            sys.exit(1)

                    if mask[mask_i] == '0':
                        new_seq += '-'
                    else:
                        new_seq += seq[seq_i]
                    mask_i += 1
                    
                seq_i += 1

            if (seq_i != seq_n):
                print("error: did not reach end of aligned sequence for taxon", tax)
                sys.exit()
            if nonmatch or (mask_i != mask_n):
                print("{0}\t{1}\t{2}\t{3}".format(tax, mask_n, mask_i, nonmatch))

            ndiff = sum([seq[i] != new_seq[i] for i in range(len(seq))])
            ndiff2 = sum([mask[i] == '0' for i in range(mask_i)])

            # ndiff = sum([seq[i] != new_seq[i] for i in range(len(seq))])
            # ndiff2 = sum([mask[i] == '0' for i in range(len(mask))])

            if ndiff != ndiff2:
                print("error in checksum on number of differences")
                print(ndiff, ndiff2)
                sys.exit()

            self.ali[tax] = new_seq

    def mask_incomplete_codons(self, noninf = "?X-*"):
        n = 0
        for (tax,seq) in self.ali.items():
            n2 = len(seq)

            if n and (n2 != n):
                print("error: sequences do not all have same length")
                sys.exit()
            else:
                n = n2

            ncodon = n2 // 3
            if 3*ncodon != n:
                print("error: sequence length not multiple of 3")
                sys.exit()
            masked_seq = ""
            for i in range(ncodon):
                if (seq[3*i] in noninf) or (seq[3*i+1] in noninf) or (seq[3*i+2] in noninf):
                    masked_seq += "---"
                else:
                    masked_seq += seq[3*i:3*(i+1)]
            self.ali[tax] = masked_seq

    def read_fasta_from_file(self, filename):
        with open(filename, 'r') as infile:
            self.read_fasta_from_stream(infile)

    def read_fasta_from_stream(self, infile):
        taxname = ''
        for line in infile:
            line = line.rstrip('\n')
            if line:
                if line[0] == '>':
                    taxname = line[1:len(line)]
                    self.ali[taxname] = ''
                else:
                    self.ali[taxname] = self.ali[taxname] + line
        self.nsite = 0
        for (tax,seq) in self.ali.items():
            if not self.nsite:
                self.nsite = len(seq)
            else:
                pass
            #    if self.nsite != len(seq):
            #        print("error in read_fasta: sequences do not have same length")
            #        sys.exit()

    def get_nonmissing_ntaxa(self, noninformative_characters = "?X-*"):
        taxlist = []
        for (tax,seq) in self.ali.items():
            informative  = sum([s not in noninformative_characters for s in seq])
            if informative:
                taxlist.append(tax)
        return len(taxlist)

    def prune_all_missing_taxa(self, noninformative_characters = "?X-*"):
        taxlist = []
        for (tax,seq) in self.ali.items():
            informative  = sum([s not in noninformative_characters for s in seq])
            if not informative:
                taxlist.append(tax)
        for tax in taxlist:
            del self.ali[tax]

    def get_missing_frac(self, noninformative_characters = "?X-*"):
        tax2miss = dict()
        for (tax,seq) in self.ali.items():
            miss = sum([s in noninformative_characters for s in seq]) / len(seq)
            tax2miss[tax] = miss
        return tax2miss

    def countmin(self, min):
        n = 0
        for (tax,seq) in self.ali.items():
            if len(seq) <= min:
                n = n + 1
        return n
        
    def prune_min_length(self, min):
        newali = dict()
        for (tax,seq) in self.ali.items():
            if len(seq) > min:
                newali[tax] = seq
        del self.ali
        self.ali = newali

    def prune_gappy_taxa(self, maxmissingfrac, noninf = "?X-*"):
        cutoff = int(maxmissingfrac * self.get_nsite())
        taxlist = []
        for (tax,seq) in self.ali.items():
            if len([s for s in seq if s not in noninf]) < min:
                taxlist.append(tax)
        for tax in taxlist:
            del self.ali[tax]

    def prune_gappy_columns(self, maxmissingfrac, noninf = "?X-*"):
        nmissing = [sum([seq[i] in noninf for tax,seq in self.ali.items()]) for i in range(self.get_nsite())]
        cutoff = int(maxmissingfrac * self.get_ntaxa())
        new_nsite = sum([m <= cutoff for m in nmissing])
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[i] for i in range(len(nmissing)) if nmissing[i] <= cutoff])
            self.ali[tax] = newseq
        self.nsite = new_nsite

    def prune_gappy_columns_per_codon(self, maxmissingfrac, noninf = "?X-*"):
        if self.get_nsite() % 3:
            print("error: not multiple of 3")
            sys.exit(1)
        ncodon = self.get_nsite() // 3
        nmissing = [sum([(seq[3*i] in noninf) or (seq[3*i+1] in noninf) or (seq[3*i+2] in noninf) for tax,seq in self.ali.items()]) for i in range(ncodon)]

        cutoff = int(maxmissingfrac * self.get_ntaxa())
        new_nsite = sum([m <= cutoff for m in nmissing])
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[3*i]+ seq[3*i+1] + seq[3*i+2] for i in range(len(nmissing)) if nmissing[i] <= cutoff])
            self.ali[tax] = newseq
        self.nsite = new_nsite

    def prune_gappy_columns_for_taxon_per_codon(self, taxname, noninf = "?X-*"):
        if self.get_nsite() % 3:
            print("error: not multiple of 3")
            sys.exit(1)
        ncodon = self.get_nsite() // 3
        seq = self.ali[taxname]
        included = [i for i in range(ncodon) if (seq[3*i] not in noninf) and (seq[3*i+1] not in noninf) and (seq[3*i+2] not in noninf)]

        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[3*i]+ seq[3*i+1] + seq[3*i+2] for i in included])
            self.ali[tax] = newseq
        self.nsite = 3*len(included)

    def prune_all_gappy_columns(self, noninf = "?X-*"):
        nmissing = [sum([seq[i] in noninf for tax,seq in self.ali.items()]) for i in range(self.get_nsite())]
        new_nsite = sum([m < self.get_ntaxa() for m in nmissing])
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[i] for i in range(len(nmissing)) if nmissing[i] < self.get_ntaxa()])
            self.ali[tax] = newseq
        self.nsite = new_nsite

    def prune_from_ali(self, other):
        toprune = [0 for i in range(self.get_nsite())]
        lastj = 0
        for i in range(other.get_nsite()):
            j = lastj
            found = 0
            while not found and j<self.get_nsite():
                if not toprune[j]:
                    diff = 0
                    for tax in self.ali:
                        diff |= self.ali[tax][j] != other.ali[tax][i]
                    if not diff:
                        found = 1
                if not found:
                    j += 1
            if j == self.get_nsite():
                print("error: could not match column from sub-alignment with complete alignment")
                sys.exit(1)
            lastj = j
            toprune[j] = 1
        if sum(toprune) != other.get_nsite():
            print("error: number of sites to prune out does not match sub-alignment length")
            sys.exit(1)
        newnsite = self.get_nsite() - other.get_nsite()
        for (tax,seq) in self.ali.items():
            newseq = "".join([seq[i] for i in range(self.get_nsite()) if not toprune[i]])
            self.ali[tax] = newseq
        self.nsite = newnsite

    def unalign(self):
        for (tax,seq) in self.ali.items():
            newseq = "".join([s for s in seq if s not in "-"])
            self.ali[tax] = newseq

    def redefine_taxon_set(self, taxon_set, all_taxa = True, noninformative_characters = "?X-*"):
        temp_tax_list = [tax for (tax,seq) in self.ali.items()]
        for tax in temp_tax_list:
            if tax not in taxon_set:
                del self.ali[tax]
        if all_taxa:
            all_missing  = "?" * self.nsite
            for tax in taxon_set:
                if tax not in self.ali:
                    self.ali[tax] = all_missing
        else:
            self.prune_all_missing_taxa(noninformative_characters = noninformative_characters)

    def taxon2genus(self):
        match_genus = r"^([A-Z][a-z]+)\_.*$"
        temp_tax_list = [tax for (tax,seq) in self.ali.items()]
        newali = dict()
        for tax in temp_tax_list:
            m = re.match(match_genus, tax)
            if not m:
                print("error in taxon2genus: could not parse taxon name {0}".format(tax))
                sys.exit(1)
            genus = m.group(1)
            if genus in newali:
                print("in taxon2genus: genus {0} appears multiple times".format(genus))
                print("over-writing previous sequence")
            newali[genus] = self.ali[tax]

        for tax in temp_tax_list:
            del self.ali[tax]

        for (genus,seq) in newali.items():
            self.ali[genus] = seq

    def write_phylip_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_phylip_to_stream(outfile)

    def write_randomized_phylip_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_randomized_phylip_to_stream(outfile)

    def write_phylip_to_stream(self, outfile):
        outfile.write("{0} {1}\n".format(self.get_ntaxa(), self.nsite))
        maxlen = 0
        for tax in sorted(self.ali.keys()):
            if maxlen < len(tax):
                maxlen = len(tax)
        for tax in sorted(self.ali.keys()):
        # for (tax,seq) in self.ali.items():
            seq = self.ali[tax]
            prettytax = tax + ''.join([' ' for i in range(maxlen - len(tax))]) 
            #outfile.write("{0}  {1}\n".format(tax,seq))
            outfile.write("{0}  {1}\n".format(prettytax,seq))

    def write_randomized_phylip_to_stream(self, outfile):
        outfile.write("{0} {1}\n".format(self.get_ntaxa(), self.nsite))
        siteranks = random.sample([i for i in range(self.get_nsite())], self.get_nsite())
        for tax in sorted(self.ali.keys()):
            seq = self.ali[tax]
            outfile.write("{0}  {1}\n".format(tax,"".join([seq[i] for i in siteranks])))

    def write_nexus_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_phylip: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_nexus_to_stream(outfile)

    def write_nexus_to_stream(self, outfile):
        outfile.write("#NEXUS\n")
        outfile.write("\n")
        outfile.write("BEGIN DATA;\n")
        outfile.write("DIMENSIONS NTAX={0} NCHAR={1};\n".format(self.get_ntaxa(), self.get_nsite()))
        outfile.write("FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
        outfile.write("MATRIX\n")
        outfile.write("\n")
        for (tax,seq) in self.ali.items():
            outfile.write("{0}  {1}\n".format(tax,seq))
        outfile.write("\t;\n")
        outfile.write("end;\n")

    def write_fasta_to_file(self, filename, force = False):
        if not force and os.path.exists(filename):
            print("in write_fasta: destination file already exists: ", filename)
            sys.exit()
        with open(filename, 'w') as outfile:
            self.write_fasta_to_stream(outfile)

    def write_fasta_to_stream(self, outfile):
        for (tax,seq) in self.ali.items():
            outfile.write(">{0}\n{1}\n".format(tax,seq))

    def sub_alignment(self, mask):
        ali = dict()
        for (tax, seq) in self.ali.items():
            newseq = ''.join([c for (i,c) in enumerate(seq) if mask[i]])
            ali[tax] = newseq
        newali = SequenceAlignment()
        newali.nsite = sum(mask)
        newali.ali = ali
        return newali

    def split(self, mask):
        ali1 = dict()
        ali2 = dict()
        for (tax, seq) in self.ali.items():
            newseq1 = ''.join([c for (i,c) in enumerate(seq) if mask[i]])
            newseq2 = ''.join([c for (i,c) in enumerate(seq) if not mask[i]])
            ali1[tax] = newseq1
            ali2[tax] = newseq2
        newali1 = SequenceAlignment()
        newali1.nsite = sum(mask)
        newali1.ali = ali1
        newali2 = SequenceAlignment()
        newali2.nsite = self.nsite - sum(mask)
        newali2.ali = ali2
        return (newali1, newali2)

    def split_equal(self, n):
        size = int(self.nsite / n)
        print(size)
        alis = []
        for k in range(n):
            ali1 = dict()
            for (tax, seq) in self.ali.items():
                newseq = seq[k*size:(k+1)*size]
                ali1[tax] = newseq
            newali = SequenceAlignment()
            newali.nsite = size
            newali.ali = ali1
            alis.append(newali)
        return alis

    def translate(self, code = univ_trans):
        if self.nsite % 3:
            print("error in translate: not multiple of 3")
            sys.exit()
        naa = self.nsite // 3
        ali = dict()
        for (tax,nucseq) in self.ali.items():
            ali[tax] = translate_nucseq(nucseq)
        protali = SequenceAlignment()
        protali.nsite = naa
        protali.ali = ali
        return protali

    def get_mean_diversity(self, noninformative_characters = "?X-*"):
        sets = [set() for i in range(self.nsite)]
        for tax,seq in self.ali.items():
            for i in range(self.nsite):
                if seq[i] not in noninformative_characters:
                    sets[i].add(seq[i])
        return sum([len(s) for s in sets])/len(sets)

    def get_gc3(self):
        meangc3 = 0
        for tax,seq in self.ali.items():
            meangc3 += get_gc3(seq)
        meangc3 /= len(self.ali)
        return meangc3

class MultiGeneSequenceAlignment:

    def __init__(self, dir_name = "", list_name = "", file_name = "", alignments = "", format = "phylip", header=False):

        self.alignments = dict()
        self.taxon_set = set()

        if list_name:
            self.read_from_list(dir_name, list_name, format = format, header= header)
        if file_name:
            self.read_from_file(file_name, format = format)
        if alignments:
            self.alignments = alignments

        self.make_taxon_set()

    def make_taxon_set(self):
        self.taxon_set = set()
        for (gene,ali) in self.alignments.items():
            self.taxon_set = self.taxon_set | ali.get_taxon_set()

    def get_taxon_set(self):
        return self.taxon_set

    def get_ntaxa(self):
        return len(self.taxon_set)

    def get_totnsite(self):
        totnsite = sum([ali.get_nsite() for gene,ali in self.alignments.items()])
        return totnsite

    def get_nsites(self):
        return [ali.get_nsite() for gene,ali in self.alignments.items()]

    def prune_min_ntaxa(self, min):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_nonmissing_ntaxa() < min:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_min_nsite(self, min):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_nsite() < min:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_max_ntaxa(self, max):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_nonmissing_ntaxa() > max:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_max_nsite(self, max):
        genelist = []
        for (gene,ali) in self.alignments.items():
            if ali.get_nsite() > max:
                genelist.append(gene)
        for gene in genelist:
            del self.alignments[gene]

    def prune_gappy_columns(self, maxmissingfrac, noninf = "?X-*"):
        for gene,ali in self.alignments.items():
            ali.prune_gappy_columns(maxmissingfrac, noninf=noninf)

    def prune_all_gappy_columns(self, noninf = "?X-*"):
        for gene,ali in self.alignments.items():
            ali.prune_all_gappy_columns(noninf=noninf)

    def prune_gappy_columns_for_taxon_per_codon(self, taxname, noninf = "?X-*"):
        for gene,ali in self.alignments.items():
            ali.prune_gappy_columns_for_taxon_per_codon(taxname, noninf=noninf)

    def get_gene_taxcount(self, min, noninf="?X-*"):
        ntaxa = self.get_ntaxa()
        gene_count = [0 for j in range(ntaxa+1)]
        for gene,ali in self.alignments.items():
            tax2miss = ali.get_missing_frac(noninformative_characters=noninf)
            ntax = ntaxa - sum([f<=min for tax,f in tax2miss.items()])
            print(gene, tax2miss, ntax)
            gene_count[ntax] += 1
        return gene_count

    def get_gene_taxcount2(self):
        ntaxa = self.get_ntaxa()
        gene_count = [0 for j in range(ntaxa+1)]
        for gene,ali in self.alignments.items():
            gene_count[ali.get_ntaxa()]+= 1
        return gene_count

    def get_mean_missing_frac(self, noninf = "?X-*"):
        tax2mean = { tax:0 for tax in self.taxon_set }

        totlength = 0
        for gene,ali in self.alignments.items():
            totlength += ali.get_nsite()
            tax2miss = ali.get_missing_frac(noninformative_characters=noninf)
            for tax in self.taxon_set:
                tax2mean[tax] += ali.get_nsite() * tax2miss[tax]
        for tax in self.taxon_set:
            tax2mean[tax] /= totlength
        return tax2mean

    def concatenate(self):
        ali = dict()
        for taxon in self.taxon_set:
            seq = "".join([self.alignments[gene].get_seq(taxon) for gene in self.alignments])
            ali[taxon] = seq
        return SequenceAlignment(fromdict=ali)

    def homogeneize_taxon_sets(self):
        for (gene,ali) in self.alignments.items():
            ali.redefine_taxon_set(self.taxon_set, all_taxa=True)

    def get_ngene(self):
        return len(self.alignments)

    def read_from_list(self, dir_name, list_name, format = "phylip", header = False):
        with open(list_name, 'r') as infile:
            ngene = 0
            if header:
                ngene = int(infile.readline().rstrip('\n'))
            gene_list = [line.rstrip('\n') for i,line in enumerate(infile) if not ngene or i<ngene]
            for gene in gene_list:
                ali = SequenceAlignment(dir_name + "/" + gene, format = format)
                self.alignments[gene] = ali

    def read_from_file(self, file_name, format = "phylip"):
        with open(file_name, 'r') as infile:
            header = infile.readline().rstrip('\n')
            if header != "ALI":
                print("error in MultiGeneSequenceAlignment.read_from_file: header")
                sys.exit()

            ngene = int(infile.readline().rstrip('\n'))

            for i in range(ngene):
                gene = ""
                while not gene:
                    gene = infile.readline().rstrip('\n')
                ali = SequenceAlignment()
                ali.read_sequential_phylip_from_stream(infile)
                self.alignments[gene] = ali
                
    def write_all_to_file(self, filename):
        with open(filename, 'w') as outfile:
            outfile.write("ALI\n")
            outfile.write("{0}\n".format(self.get_ngene()))
            for (gene,ali) in self.alignments.items():
                outfile.write("{0}\n".format(gene))
                ali.write_phylip_to_stream(outfile)

    def write_all_to_files(self, basename, force = False):
        for gene,ali in self.alignments.items():
            filename = basename + gene
            if not force and os.path.exists(filename):
                print("error when writing gene-specific alignments: files {0} already exists".format(filename))
                sys.exit(1)
            with open(filename, 'w') as outfile:
                ali.write_phylip_to_stream(outfile)

    def write_genes_to_nexus(self, dirname):
        for (gene,ali) in self.alignments.items():
            ali.write_nexus_to_file(dirname + "/" + gene)

    def get_gene(self, gene):
        return self.alignments[gene]

    def get_gene_list(self):
        gene_list = sorted([gene for (gene,ali) in self.alignments.items()])
        return gene_list

    def get_gene_ali(self, gene):
        if gene not in self.alignments:
            print("error: gene not in multigene sequence alignment")
            sys.exit(1)
        return self.alignments[gene]

    def change_gene_names(self, name_pattern):
        alignments = dict()
        for (gene, ali) in self.alignments:
            m = re.match(pattern, gene)
            if not m:
                print("error in MultiGeneSequenceAlignment.read_from_list: gene name does not match pattern")
                sys.exit()
            new_name = m.group(1)
            alignments[new_name] = self.alignments[gene]
        self.alignments = alignments

    def prune_all_missing_taxa(self, noninformative_characters = "?X-*"):
        for (gene,ali) in self.alignments.items():
            ali.prune_all_missing_taxa(noninformative_characters = noninformative_characters)

    def redefine_taxon_set(self, taxon_set, all_taxa = True, noninformative_characters = "?X-*"):
        for (gene,ali) in self.alignments.items():
            ali.redefine_taxon_set(taxon_set, all_taxa = all_taxa, noninformative_characters = noninformative_characters)
        self.make_taxon_set()

    def translate(self, code = univ_trans):
        alignments = dict()
        for (gene,ali) in self.alignments.items():
            protali = ali.translate(code = code)
            alignments[gene] = protali
        prot_multiali = MultiGeneSequenceAlignment(alignments=alignments)
        return prot_multiali

    def get_mean_diversity(self, gene, noninformative_characters="?X-*"):
        return self.alignments[gene].get_mean_diversity(noninformative_characters = noninformative_characters)

    def get_gc3(self, gene):
        return self.alignments[gene].get_gc3()

