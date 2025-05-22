#!/usr/bin/env python3
import re
import argparse

def find_polyA_from_left(seq, min_len=12):
    """
    Find the longest A stretch from the left (5' to 3') of the sequence.

    Returns:
        start (int): start index of polyA
        end (int): end index (exclusive)
        polya_seq (str): the polyA sequence
    """
    match = re.match(r'(A{' + str(min_len) + r',})', seq)
    if match:
        start = match.start()
        end = match.end()
        polya_seq = match.group()
        return start, end, polya_seq
    else:
        return None, None, None


def main():
    parser = argparse.ArgumentParser(description="Find longest polyA stretch from the left side of the sequence")
    parser.add_argument('--seq', type=str, required=True, help='Input sequence to scan for polyA')
    parser.add_argument('--min_len', type=int, default=12, help='Minimum polyA length (default: 12)')
    
    args = parser.parse_args()
    seq = args.seq.upper()

    start, end, polya = find_polyA_from_left(seq, args.min_len)
    
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
        print(f"No polyA stretch of at least {args.min_len} bases found at the left of the sequence.")

if __name__ == "__main__":
    main()
