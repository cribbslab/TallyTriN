#!/usr/bin/env python3

import argparse
import parasail
import gzip
import pysam

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

def process_sequence(seq):
    seq_tail = seq[-200:]
    result = is_adapter_match(seq_tail)
    print(f"\n--- Result ---")
    print(f"Tail (200bp): {result['tail']}")
    print(f"Adapter: {ADAPTER}")
    print(f"Adapter Match: {result['matched_seq']}")
    print(f"Aligned Seq: {result['aligned_seq']}")
    print(f"Indices: {result['start_idx']} to {result['end_idx']}")
    print(f"Score: {result['score']:.2f} (Threshold: {MIN_SCORE:.2f})")
    print(f"Match: {'PASS' if result['match'] else 'FAIL'}")
    print(f"UMI Window: {result['umi_window'] if result['umi_window'] else 'N/A'}")

def main():
    parser = argparse.ArgumentParser(description="Test adapter match in sequence or FASTQ.")
    parser.add_argument("--seq", help="Sequence string or FASTQ/FASTQ.GZ file path")
    parser.add_argument("-n", "--num", type=int, default=1, help="Number of reads to process")
    args = parser.parse_args()

    if args.seq.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        with pysam.FastxFile(args.seq) as fh:
            for i, record in enumerate(fh):
                if i >= args.num:
                    break
                print(f"\n== Read {i+1}: {record.name} ==")
                process_sequence(record.sequence)
    else:
        # Assume raw sequence string
        process_sequence(args.seq)

if __name__ == "__main__":
    main()
