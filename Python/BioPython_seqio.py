#!/usr/bin/env python3
"""Question 4"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# pylint: disable=invalid-name


def write_to_file(infile, outfile):
    """Takes infile, prints out reverse
    complement to outfile."""
    # open outfile so that it closes on its own
    with open(outfile, 'w') as outfile:
        for seq_record in SeqIO.parse(infile, 'fasta'):
            my_seq = SeqRecord(id=seq_record.id,
                               seq=seq_record.seq.reverse_complement())
            SeqIO.write(my_seq, outfile, 'fasta')


if __name__ == "__main__":

    # check number of arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: an infile and "
                        "an outfile")

    # takes two arguments from command line: original FASTA and new FASTA

    infile = sys.argv[1]
    outfile = sys.argv[2]

    write_to_file(infile, outfile)
