import numpy as np
import sys
import os

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Align import write
from Bio.SeqIO.FastaIO import SimpleFastaParser

def main():

    orig_fasta = sys.argv[1]
    dir_name = os.path.dirname(sys.argv[1])

    seqs_list = []
    with open(orig_fasta, 'r') as handle:
        for values in SimpleFastaParser(handle):
            seqs_list.append(values[1])

    pairs = []
    for i in range(len(seqs_list)):
        for j in range(i+1, len(seqs_list)):
            pairs.append((seqs_list[i], seqs_list[j]))

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = sys.argv[2]
    aligner.gap_score = -0.5

    for index, pair in enumerate(pairs):
        alignments = aligner.align(pair[0], pair[1])
        #print(aligner.algorithm)
        alignment = next(alignments)
        mode = sys.argv[2]
        file_name = dir_name + '/' + mode + '-' + str(index) + '.fasta'
        with open(file_name, 'w') as handle:
            handle.write(alignment.format("fasta"))
    

if __name__ == "__main__":
    main()