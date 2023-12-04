import os
import sys

import numpy as np
import pandas as pd

from Bio.SeqIO.FastaIO import SimpleFastaParser

BLOSUM_MATRIX_PATH = 'blosum62.csv'

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
               'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z']


def generate_dict_conservative_mutation():
    blosum = pd.read_csv(BLOSUM_MATRIX_PATH)
    blosum = blosum.to_numpy()
    dict_conservative = {}
    i = 0
    for row in blosum:
        positives = np.where(row > 0)
        positives = positives[0].tolist()
        positives_aa = [amino_acids[pos] for pos in positives]
        dict_conservative[amino_acids[i]] = positives_aa
        i += 1

    return dict_conservative


def generate_alignment(seqs_list, conserv_dict):
    seq1 = seqs_list[0]
    seq2 = seqs_list[1]
    identities = 0
    conservations = 0
    gaps = 0
    alignment = ''
    for index in range(len(seq1)):
        if seq1[index] == '-' or seq2[index] == '-':
            gaps += 1
            alignment += '-'
        elif seq1[index] == seq2[index]:
            identities += 1
            conservations += 1
            alignment += seq1[index]
        elif seq1[index] in conserv_dict[seq2[index]]:
            conservations += 1
            alignment += '+'
        else:
            alignment = ' '

    return (identities/len(seq1), conservations/len(seq1), gaps/len(seq1), alignment)


def get_sequences_list(file_name):
    seqs_list = []
    dir_file = sys.argv[1] + '/' + file_name
    with open(dir_file, 'r') as handle:
        for values in SimpleFastaParser(handle):
            seqs_list.append(values[1])

    return seqs_list


def main():

    global_aligns = []
    local_aligns = []

    files = os.listdir(sys.argv[1])
    dict_cons = generate_dict_conservative_mutation()

    for fil in files:
        if fil[0:6] == 'global':
            global_aligns.append(fil)
        elif fil[0:5] == 'local':
            local_aligns.append(fil)
        else:
            pass

    # Global evaluation
    id_glob = []
    cons_glob = []
    gap_glob = []

    for glob in global_aligns:
        seqs_list = get_sequences_list(glob)
        (id, cons, gap, _) = generate_alignment(seqs_list, dict_cons)
        id_glob.append(id)
        cons_glob.append(cons)
        gap_glob.append(gap)

    id_glob = np.array(id_glob)
    cons_glob = np.array(cons_glob)
    gap_glob = np.array(gap_glob)

    print('Global Alignment Average Scores:')
    print(
        f'Identity: {np.mean(id_glob)} - Conservation: {np.mean(cons_glob)} - Gaps: {np.mean(gap_glob)}')
    print('----------------------------------------------------------------------------------------')

    # Local evaluation
    id_loc = []
    cons_loc = []
    gap_loc = []

    for loc in local_aligns:
        seqs_list = get_sequences_list(loc)
        (id, cons, gap, _) = generate_alignment(seqs_list, dict_cons)
        id_loc.append(id)
        cons_loc.append(cons)
        gap_loc.append(gap)

    id_loc = np.array(id_loc)
    cons_loc = np.array(cons_loc)
    gap_loc = np.array(gap_loc)

    print('Local Alignment Average Scores:')
    print(
        f'Identity: {np.mean(id_loc)} - Conservation: {np.mean(cons_loc)} - Gaps: {np.mean(gap_loc)}')
    print('----------------------------------------------------------------------------------------')


if __name__ == '__main__':
    main()
