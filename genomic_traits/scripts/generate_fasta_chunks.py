#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
import os

def batch_iterator(iterator, chunk_size):
    """Returns lists of length chunk_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    chunk_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == chunk_size:
            yield batch
            batch = []
    if batch:
        yield batch


def create_batch_fasta(fpath, odpath, chunk_size):
    basefname, extension = os.path.splitext(os.path.basename(relative_path))


    
    record_iter = SeqIO.parse(open(fpath), "fasta")
    for i, batch in enumerate(batch_iterator(record_iter, chunk_size), start=1):
        out_fpath  = os.path.join(odpath, f"{basefname}_{i}{extension}") 
        with open(out_fpath, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, out_fpath))



if __name__ == '__main__':
    import argparse
    import json
    import pprint

    parser = argparse.ArgumentParser(description='Run models - nutrients recycle with separate N/C and quotas.')
    parser.add_argument('--fasta', help='fasta to split', required=True)
    parser.add_argument("--out_dpath", help="output dir", default='.')
    parser.add_argument("--chunk_size", help="chunk_size",  type=int, default=100)
          
    
    args = parser.parse_args()
    odpath = args.out_dpath
    if odpath != '':
        os.makedirs(odpath, exist_ok=True)

    create_batch_fasta(args.fasta, odpath, args.chunk_size)
    