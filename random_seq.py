import argparse
import random

import numpy as np
import pandas as pd

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
               'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z']

amino_acids_dict = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7,
                    'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                    'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20,
                    'Z': 21}

MAX_MUTATION_RATE = 0.6


def create_random_sequence(seq_size):
    random_seq = [random.choice(amino_acids) for _ in range(seq_size)]
    return ''.join(random_seq)


def create_blosum_matrix(csv):
    blosum = pd.read_csv(csv)
    blosum = blosum.to_numpy()

    return blosum


def mutate_sequence(original_seq, mutation_rate, blosum, epsilon, indel, indel_rate):
    new_seq = list(original_seq)

    num_mutations = int(len(original_seq) * mutation_rate)
    mutated_indexes = random.sample(range(len(original_seq)), num_mutations)

    for index in mutated_indexes:
        amino_acid = original_seq[index]
        mutated_amino_acid = get_blosum_choice(amino_acid, blosum, epsilon)
        new_seq[index] = mutated_amino_acid

    if indel:
        num_indels = int(len(original_seq) * indel_rate)
        indel_indexes = random.sample(range(len(original_seq)), num_indels)

        for index in indel_indexes:
            if index >= len(new_seq):
                index = len(new_seq) - 1

            if random.random() < 0.5:
                # Insertion
                insertion = random.choice(amino_acids)
                new_seq.insert(index, insertion)
            else:
                # Deletion
                new_seq.pop(index)

    return ''.join(new_seq)


def get_blosum_choice(amino_acid, blosum, epsilon):
    amino_acid_index = amino_acids_dict[amino_acid]
    sorted_scores = np.argsort(blosum[amino_acid_index])

    if random.random() < epsilon:
        index = random.choice(sorted_scores[-4:-1])
    else:
        index = random.choice(sorted_scores[:-4])

    return amino_acids[index]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_size', type=int,
                        help='Size of the original random sequence')
    parser.add_argument('num_mutated_seqs', type=int,
                        help='Number of mutated sequences to be produced')
    parser.add_argument('mutation_rate', type=float,
                        help='Rate of occurrence of mutations')
    parser.add_argument('epsilon', type=float,
                        help='Epsilon probability used for mutation')
    parser.add_argument('--indel', action='store_true', default=False,
                        help='Determines if mutations should have indels')
    parser.add_argument('--progressive', action='store_true', default=False,
                        help='Determines if mutation rates should increase')
    parser.add_argument('--indel_total_rate', type=float,
                        default=0.1, help='Rate of sequences that should have indels')
    parser.add_argument('--indel_rate', type=float,
                        default=0.01, help='Rate of indels with respect to the size of the sequence')
    args = parser.parse_args()

    blosum_matrix = create_blosum_matrix('blosum62.csv')
    original_seq = create_random_sequence(args.seq_size)

    print('>Original')
    print(f'{original_seq}')

    if args.indel:
        num_indel_sequences = int(
            args.num_mutated_seqs * args.indel_total_rate)
        indel_seqs = set(random.sample(
            range(args.num_mutated_seqs), num_indel_sequences))
    else:
        indel_seqs = set()

    mutants = []
    if not args.progressive:
        for index in range(args.num_mutated_seqs):
            is_indel = index in indel_seqs
            mutant = mutate_sequence(original_seq, args.mutation_rate,
                                blosum_matrix, args.epsilon, is_indel, args.indel_rate)
            mutants.append(mutant)
    else:
        rates = np.linspace(args.mutation_rate,
                            MAX_MUTATION_RATE, args.num_mutated_seqs)
        for (index, mutation_rate) in enumerate(rates):
            is_indel = index in indel_seqs
            mutant = mutate_sequence(original_seq, mutation_rate,
                                blosum_matrix, args.epsilon, is_indel, args.indel_rate)
            mutants.append(mutant)
        
    random.shuffle(mutants)
    for (i, mutant) in enumerate(mutants):
        print(f'>Mutant{i}')
        print(mutant)


if __name__ == '__main__':
    main()
