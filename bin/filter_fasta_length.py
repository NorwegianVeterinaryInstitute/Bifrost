#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def main(input, output, minimum):
    with open(input, "rU") as input_handle, open(output, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record) >= minimum:
                SeqIO.write(record, output_handle, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", metavar="STRING", help = "input filename")
    parser.add_argument("-o", "--output", metavar="STRING", help = "output filename")
    parser.add_argument("-m", "--minimum", metavar="INT", help = "minimum contig length \
                                                        (smaller than will be removed)")
    args = parser.parse_args()

    main(args.input, args. output, int(args.minimum))
