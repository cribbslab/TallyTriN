#!/usr/bin/env python3

import argparse
import parasail
import pysam
import os

# Constants
ADAPTER = "GTACTCTGCGTTGATACCACTGCTT"
MATCH_SCORE = 2
MISMATCH_PENALTY = -1
ERROR_RATE = 0.35
PRIMER_LENGTH = len(ADAPTER)
MIN_SCORE = (1 - ERROR_RATE) * PRIMER_LENGTH * MATCH_SCORE
MATRIX = parasail.matrix_create("ACGT", MATCH_SCORE, MISMATCH_PENALTY)

def is_adapter_match(seq_tail):
    result = parasail.sw_trace_striped_16(ADAPTER, seq_tail, 1, 1, MATRIX)
    score = result.score
    aligned_len = len(result.traceback.ref.replace("-", ""))
    start_idx = result.end_ref - aligned_len + 1
    end_idx = result.end_ref + 1
    matched_seq = seq_tail[start_idx:end_idx]
    aligned_seq = result.traceback.ref
    match_pass = score >= MIN_SCORE
    umi_start = start_idx - 60
    umi_window = seq_tail[umi_start:start_idx] if umi_start >= 0 else None
    return {
        "tail": seq_tail,
        "matched_seq": matched_seq,
        "start_idx": start_idx,
        "end_idx": end_idx,
        "aligned_seq": aligned_seq,
        "score": score,
        "match": match_pass,
        "umi_window": umi_window
    }

def process_sequence(seq, label="sequence"):
    seq_tail = seq[-200:]
    result = is_adapter_match(seq_tail)
    print(f"\n== {label} ==")
    print(f"Tail (200bp): {result['tail']}")
    print(f"Adapter: {ADAPTER}")
    print(f"Adapter Match: {result['matched_seq']}")
    print(f"Aligned Seq: {result['aligned_seq']}")
    print(f"Indices: {result['start_idx']} to {result['end_idx']}")
    print(f"Score: {result['score']:.2f} (Threshold: {MIN_SCORE:.2f})")
    print(f"Match: {'PASS' if result['match'] else 'FAIL'}")
    print(f"UMI Window: {result['umi_window'] if result['umi_window'] else 'N/A'}")

def iter_records_from_file(path, n=1):
    with pysam.FastxFile(path) as fh:
        for i, record in enumerate(fh):
            if i >= n:
                break
            yield record.name, record.sequence

def main():
    parser = argparse.ArgumentParser(description="Test adapter match in a raw sequence or FASTA/FASTQ file.")
    parser.add_argument("--seq", help="Raw sequence string")
    parser.add_argument("--infile", help="Input .fa/.fasta/.fq/.fastq file (optionally gzipped)")
    parser.add_argument("-n", "--num", type=int, default=1, help="Number of reads to process from file input")
    args = parser.parse_args()

    if args.seq and args.infile:
        raise ValueError("Use either --seq or --infile, not both.")
    if not args.seq and not args.infile:
        raise ValueError("You must provide either --seq or --infile.")

    if args.infile:
        for i, (name, seq) in enumerate(iter_records_from_file(args.infile, args.num), 1):
            process_sequence(seq, label=f"read {i}: {name}")
    else:
        process_sequence(args.seq, label="raw sequence")

if __name__ == "__main__":
    main()
