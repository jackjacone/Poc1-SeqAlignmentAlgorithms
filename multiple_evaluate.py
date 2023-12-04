from pair_evaluate import generate_dict_conservative_mutation

import argparse
import itertools

from Bio import AlignIO
import numpy as np


def compute_metrics(alignment):
    conservation_dict = generate_dict_conservative_mutation()
    identities = 0
    gaps = 0
    conservations = 0
    length = alignment.get_alignment_length()

    for idx in range(length):
        s = alignment[:, idx]

        if '-' in s:
            gaps += 1
            continue

        if s == len(s) * s[0]:
            identities += 1
            conservations += 1
        else:
            unique = set(s)
            if all(a1 in conservation_dict[a2] for (a1, a2) in itertools.combinations(unique, 2)):
                conservations += 1

    identity = identities/length
    conservation = conservations/length
    gaps = gaps/length

    return identity, conservation, gaps


def synthetic_conservation(alignment):
    conservation_dict = generate_dict_conservative_mutation()

    original_seq = alignment[0]
    conservations_list = []
    for i in range(1, len(alignment)):
        conservations = 0
        for (j, value) in enumerate(alignment[i]):
            original_value = original_seq[j]
            if original_value == '-' or value == '-':
                continue

            if value in conservation_dict[original_value]:
                conservations += 1

        conservations_list.append(
            conservations/alignment.get_alignment_length())

    return np.mean(conservations_list)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment_fasta', type=str,
                        help='Path to the alignment\'s FASTA file')
    parser.add_argument('--is_synthetic', action='store_true', default=False,
                        help='Determines if the alignment is synthetic')
    args = parser.parse_args()

    with open(args.alignment_fasta, 'r') as fasta_file:
        alignment = AlignIO.read(fasta_file, 'fasta')
        identity, conservation, gaps = compute_metrics(alignment)

        if args.is_synthetic:
            conservation = synthetic_conservation(alignment)

        print(
            f'Identity: {identity} - Conservation: {conservation} - Gaps: {gaps}')


if __name__ == '__main__':
    main()
