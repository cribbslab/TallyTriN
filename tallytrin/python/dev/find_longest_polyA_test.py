#!/usr/bin/env python3
import re
import argparse
import os
import gzip

def find_polyA_from_left(seq, min_len=12):
    """
    Find the longest A stretch from the left (5' to 3') of the sequence.
    """
    match = re.search(r'(A{' + str(min_len) + r',})', seq)
    if match:
        start = match.start()
        end = match.end()
        polya_seq = match.group()
        return start, end, polya_seq
    else:
        return None, None, None

def read_sequences(path, n=1):
    """
    Yield up to n sequences from FASTA/FASTQ (plain or gzipped).
    For FASTQ, sequence is line 2 of each 4-line record.
    For FASTA, concatenates sequence lines until next header.
    """
    opener = gzip.open if path.endswith(".gz") else open
    ext = path.lower().replace(".gz", "")

    with opener(path, "rt") as fh:
        if ext.endswith((".fq", ".fastq")):
            while n != 0:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()  # +
                fh.readline()  # quality
                yield header.strip()[1:] if header.startswith("@") else "read", seq
                n -= 1
        elif ext.endswith(".fa") or ext.endswith(".fasta"):
            name = None
            seq_chunks = []
            count = 0
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(seq_chunks)
                        count += 1
                        if n != 0 and count >= n:
                            return
                    name = line[1:]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if name is not None and (n == 0 or count < n):
                yield name, "".join(seq_chunks)
        else:
            raise ValueError("Unsupported file type. Use .fa/.fasta/.fq/.fastq, optionally .gz")

def process_sequence(seq, label="sequence", min_len=12):
    seq = seq.upper()
    start, end, polya = find_polyA_from_left(seq, min_len)

    print(f"\n=== {label} ===")

    if polya:
        print(f"PolyA stretch found: {polya}")
        print(f"Start position: {start}")
        print(f"End position: {end}")
        print(f"Length of polyA: {len(polya)}")
        print(f"Base after polyA (boundary base): {seq[end] if end < len(seq) else 'END'}")
        print(f"Remaining sequence after polyA: {seq[end+1:]}")

        umi_seq = seq[end+1:]

        if umi_seq:
            print(f"UMI Candidate Length: {len(umi_seq)}")
        else:
            print("\n❌ Could not find valid UMI.")
    else:
        print(f"No polyA stretch of at least {min_len} bases found at the left of the sequence.")

def main():
    parser = argparse.ArgumentParser(
        description="Find longest polyA stretch from the left side of a sequence or file"
    )
    parser.add_argument('--seq', type=str, default=None,
                        help='Raw input sequence')
    parser.add_argument('--infile', type=str, default=None,
                        help='Input FASTA/FASTQ file (.fa/.fq/.fastq, optionally .gz)')
    parser.add_argument('--n', type=int, default=1,
                        help='Number of reads to process from file input (default: 1)')
    parser.add_argument('--min_len', type=int, default=12,
                        help='Minimum polyA length (default: 12)')

    args = parser.parse_args()

    if args.seq and args.infile:
        raise ValueError("Use either --seq or --infile, not both.")
    if not args.seq and not args.infile:
        raise ValueError("You must provide either --seq or --infile.")

    if args.seq:
        process_sequence(args.seq, label="input_sequence", min_len=args.min_len)
    else:
        for i, (name, seq) in enumerate(read_sequences(args.infile, args.n), 1):
            process_sequence(seq, label=f"{name} (read {i})", min_len=args.min_len)

if __name__ == "__main__":
    main()
