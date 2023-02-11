#!/usr/bin/env python3
"""Question 5"""

import sys
from Bio import SeqIO

# pylint: disable=redefined-outer-name
# pylint: disable=consider-using-enumerate
# pylint: disable=invalid-name


def sliding_window(k, string):
    """ takes k-mer size and string, returns list
    of all kmers in the given string using the
    sliding window algorithm.
    :param k:
    :param string:
    :return: kmers
    """
    # initiate kmers list
    kmers = []

    end = len(string) - k + 1
    for start in range(0, end):
        kmers.append(string[start:start + k])

    return kmers


def gc_content(dna):
    """Takes single string and returns it's GC content
    as a fraction between 0 and 1.
    :param dna: string
    :return:  gc_1 / len(dna)
    """
    # change dna to lowercase to make it more universal
    dna = dna.lower()
    gc_1 = 0

    for nucleotide in dna:
        if nucleotide in ['g', 'c']:
            gc_1 += 1

    return gc_1 / len(dna)


if __name__ == "__main__":

    # check number of arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: a kmer size and "
                        "then a fasta file")

    # define arguments
    k = int(sys.argv[1])
    fasta_file = sys.argv[2]

    # parse fasta file
    my_seq = SeqIO.parse(fasta_file, 'fasta')

    # print out description and sequences
    for item in my_seq:
        print(f">{item.description}")
        for kmer in sliding_window(k, item.seq):
            print("{}\t{:9.2f}".format(kmer, gc_content(kmer)))
