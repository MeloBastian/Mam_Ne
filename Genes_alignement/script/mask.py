import sys
from ali import *

aliname = sys.argv[1]
maskname = sys.argv[2]
unaliname = sys.argv[3]
outname = sys.argv[4]

ali = SequenceAlignment(aliname, format = 'fasta')
unali = SequenceAlignment(unaliname, format = 'fasta')
mask = SequenceAlignment(maskname, format = 'fasta')
ali.mask_ali(mask, unali)
ali.mask_incomplete_codons()
ali.write_fasta_to_file(outname)

