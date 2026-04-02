import sys
import regex
import cgatcore.iotools as iotools
import pysam
import logging
import argparse
import collections
from collections import Counter
import parasail

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("polya_umi.py")

# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--infile", default=None, type=str,
                    help='infile fastq  file')
parser.add_argument("--outname", default=None, type=str,
                    help='name for output fastq files')
parser.add_argument("--errors", default=None, type=str,
                    help='Number of errors to remove reads')
args = parser.parse_args()

L.info("args:")
print(args)



# ########################################################################### #
# ######################## Code                ############################## #
# ########################################################################### #

# ------------------ Helper Functions ------------------

def find_longest_polyA(seq):
    max_len = 0
    max_start = 0
    current_len = 0
    current_start = 0

    for i, base in enumerate(seq):
        if base == 'A':
            if current_len == 0:
                current_start = i
            current_len += 1
            if current_len > max_len:
                max_len = current_len
                max_start = current_start
        else:
            current_len = 0

    max_end = max_start + max_len - 1
    return max_start, max_end, max_len

def extract_umi_from_polyA_region(seq):
    start, end, polyA_len = find_longest_polyA(seq)

    if polyA_len < 12:
        return None, 0, "", 0, 0, 0

    boundary_index = end + 2
    if boundary_index >= len(seq):
        return None, 0, "", 0, 0, 0

    remaining = seq[boundary_index:]
    if len(remaining) < 30 or len(remaining) > 47:
        return None, 0, remaining, 0, 0, 0
      
    max_errors = {'ins': 10, 'del': 10, 'subs': 10}
    max_total_errors = 10

    def is_homotrimer(trimer):
        return len(trimer) == 3 and trimer[0] == trimer[1] == trimer[2]
    
    def is_soft_homotrimer(trimer):
        # Returns True if the trimer has at least 2 identical bases
        if len(trimer) != 3:
            return False
        base_counts = Counter(trimer)
        return base_counts.most_common(1)[0][1] >= 2

    def get_disrupting_position(trimer, target_base):
        for i, b in enumerate(trimer):
            if b != target_base:
                return i
        return None

 
    def try_deletion_strategies(seq, trimer_start, extra_bases):
        deletions = 0
        attempts = []
    
        # Prepare positions for deletions
        positions = [[0], [1], [2], [0, 1], [1, 2], [2, 3]]  # 3 single, 3 double deletions
    
        for pos_list in positions:
            temp_seq = seq
            del_indices = [trimer_start + p for p in pos_list if trimer_start + p < len(temp_seq)]
    
            if len(del_indices) != len(pos_list):
                continue  # Skip if any deletion index is out of range
    
            removed_bases = ''.join(temp_seq[i] for i in del_indices)
            # Perform deletion
            for offset, i in enumerate(del_indices):
                adjusted_i = i - offset  # adjust index as sequence shrinks
                temp_seq = temp_seq[:adjusted_i] + temp_seq[adjusted_i + 1:]
    
            # Extract 9bp window after deletion
            window_start = trimer_start
            temp_window = temp_seq[window_start:window_start + 9]
    
            # Pad with extra bases if available
            needed = 9 - len(temp_window)
            if needed > 0 and len(extra_bases) >= 0:
                temp_window += extra_bases[:needed]
    
            # Break into trimers and test
            trimers = [temp_window[i:i + 3] for i in range(0, len(temp_window), 3)]
    
            #print(f"\n🧪 Deletion attempt: positions {pos_list} (removed: {removed_bases})")
            #print(f"🔄 New 9bp window: {temp_window}")
            #print(f"🧬 Trimers: {trimers}")
            
            if len(trimers) > 0 and all(len(t) == 3 and is_homotrimer(t) for t in trimers):
            #if len(trimers) == 3 and all(len(t) == 3 and t[0] == t[1] == t[2] for t in trimers):
                # If successful, patch this corrected window into the full sequence
                corrected_seq = temp_seq
    
                # Reconstruct full 36bp window using extra_bases
                if len(corrected_seq) < UMI_LENGTH:
                    needed_full = UMI_LENGTH - len(corrected_seq)
                    if len(extra_bases) >= needed_full:
                        corrected_seq += extra_bases[:needed_full]
                    else:
                        #print("❌ Not enough extra bases to rebuild full 36bp window after deletion.")
                        continue
    
                corrected_seq = corrected_seq[:UMI_LENGTH]  # ensure exact length
    
                #print("✅ Deletion successful.")
                return corrected_seq, len(pos_list), 'del'
    
        #print("❌ All deletion attempts failed.")
        return None, 0, ''

 
    def try_insertion_strategies(seq, trimer_start):
        trimer = seq[trimer_start:trimer_start + 3]
        base_counts = Counter(trimer)
        most_common = base_counts.most_common(1)[0][0]
        least_common = base_counts.most_common()[-1][0]

        attempts = [
            (most_common, 1),
            (most_common, 2),
            (least_common, 1),
            (least_common, 2)
        ]

        for base, count in attempts:
            pos = get_disrupting_position(trimer, most_common)
            if pos is None:
                pos = 1
            insert_pos = trimer_start + pos
            new_seq = seq[:insert_pos] + base * count + seq[insert_pos:]
            window = new_seq[trimer_start:trimer_start + 9]
            trimers = [window[i:i + 3] for i in range(0, len(window), 3)]
            #print(f"Inserting '{base * count}' at {insert_pos}, new window: {window}")
            #print(f"Trimers after insertion: {trimers}")
            if len(trimers) > 0 and all(len(t) == 3 and is_homotrimer(t) for t in trimers):
                #print("Insertion successful.")
                return new_seq, count, 'ins'
        #print("All insertion attempts failed.")
        return None, 0, ''
      

    def try_substitution_strategies(modified_seq, i, trimer, strict, used_errors_bases, used_errors_events):
        # Strict mode requires next two trimers to be perfect homotrimers
        if strict:
            next1 = modified_seq[i+3:i+6]
            next2 = modified_seq[i+6:i+9]
            if not (is_homotrimer(next1) and is_homotrimer(next2)):
                return None, 0, ''
        else:
            # Relaxed mode: Check fallback window for "soft homotrimers"
            fallback_window = modified_seq[i:i+9]
            fallback_trimers = [fallback_window[j:j+3] for j in range(0, len(fallback_window) - 2, 3)]
            if not (len(fallback_trimers) > 0 and all(is_soft_homotrimer(t) for t in fallback_trimers)):
               return None, 0, ''
    
        base_counts = Counter(trimer)
        if len(base_counts) == 2 and max(base_counts.values()) == 2:
            corrected_base = max(base_counts, key=base_counts.get)
            corrected_trimer = corrected_base * 3
            modified_seq = modified_seq[:i] + corrected_trimer + modified_seq[i+3:]
            return modified_seq, 1, 'subs'
    
        return None, 0, ''


    def correct_and_validate(seq, extra_bases, max_errors=None, max_total_errors=None):
        def threshold_exceeded():
            return (max_total_errors is not None and total_errors > max_total_errors) or \
                   (max_errors and any(used_errors_bases[k] > max_errors[k] for k in max_errors))
        
        
        def track_error(err_type, errors=1):
            nonlocal total_errors, corrected_trimers
            total_errors += errors
            corrected_trimers += 1
            used_errors_bases[err_type] += errors
            used_errors_events[err_type] += 1
        
        modified_seq = seq
        i = 0
        total_errors = 0
        corrected_trimers = 0  # track how many trimers needed correction
        used_errors_bases = {'ins': 0, 'del': 0, 'subs': 0}
        used_errors_events = {'ins': 0, 'del': 0, 'subs': 0}
    
        while i < len(modified_seq):
            trimer = modified_seq[i:i + 3]
            if is_homotrimer(trimer):
                i += 3
                continue
            # === Final trimer at position i attempting to fix via insertion ===
            if len(trimer) < 3:
                new_seq, errors, err_type = try_insertion_strategies(modified_seq, i)
                if new_seq:
                    modified_seq = new_seq
                    track_error(err_type, errors)
                    # Threshold check
                    if threshold_exceeded():
                        return None, float('inf'), 0, used_errors_bases, used_errors_events
                    continue
                  
                return None, float('inf'), corrected_trimers, used_errors_bases, used_errors_events
            
            # === Substitution Check (strict): needs 3x 2/3 match in 9bp window ===
            new_seq, errors, err_type = try_substitution_strategies(
                modified_seq, i, trimer, strict=True, 
                used_errors_bases=used_errors_bases, 
                used_errors_events=used_errors_events)
    
            if new_seq:
                modified_seq = new_seq
                track_error(err_type, errors)
                
                i += 3
                # Threshold check
                if threshold_exceeded():
                    return None, float('inf'), 0, used_errors_bases, used_errors_events
                continue
              
            # === Try indel strategies ===    
            tried_del = tried_ins = False
            if len(modified_seq) >= UMI_LENGTH:
                new_seq, errors, err_type = try_deletion_strategies(modified_seq, i, extra_bases)
                if new_seq:
                    tried_del = True
            else:
                new_seq, errors, err_type = try_insertion_strategies(modified_seq, i)
                if new_seq:
                    tried_ins = True

           
            if not tried_del and not tried_ins:
                if len(modified_seq) >= UMI_LENGTH:
                    new_seq, errors, err_type = try_insertion_strategies(modified_seq, i)
                    if new_seq:
                        tried_ins = True
                else:
                    new_seq, errors, err_type = try_deletion_strategies(modified_seq, i, extra_bases)
                    if new_seq:
                        tried_del = True
    
            if tried_del or tried_ins:
                modified_seq = new_seq
                track_error(err_type, errors)
                # Threshold check
                if threshold_exceeded():
                    return None, float('inf'), 0, used_errors_bases, used_errors_events
                continue
 
           # === Fallback substitution (if all indels failed and all 3 have 2/3 match) ===

            new_seq, errors, err_type = try_substitution_strategies(
                modified_seq, i, trimer, strict=False,
                used_errors_bases=used_errors_bases,
                used_errors_events=used_errors_events)
            if new_seq:
                modified_seq = new_seq
                track_error(err_type, errors)
                i += 3
                # Threshold check
                if threshold_exceeded():
                    #print(f"❌ Threshold exceeded: {used_errors_bases}, abandoning window.")
                    return None, float('inf'), 0, used_errors_bases, used_errors_events
                continue


            return None, float('inf'), 0, used_errors_bases, used_errors_events

        # === Final validation step: 12 perfect homotrimers ===    
        trimers = [modified_seq[i:i + 3] for i in range(0, len(modified_seq), 3)]
    
        if len(trimers) >= UMI_TRIMERS and all(is_homotrimer(t) for t in trimers[:UMI_TRIMERS]):
            umi = "".join(trimers[:UMI_TRIMERS])
            excess = len(umi) - UMI_LENGTH
            input_overhang = len(seq) - UMI_LENGTH if len(seq) > UMI_LENGTH else 0
            total_excess = max(excess, input_overhang)
    
            if total_excess > 0:
               total_errors += total_excess
               umi = umi[:UMI_LENGTH]
    
            return umi, total_errors, corrected_trimers, used_errors_bases, used_errors_events
        else:
            return None, float('inf'), corrected_trimers, used_errors_bases, used_errors_events
      
    
    # ==== MAIN BODY ====
    #print(f"Initial sequence: {seq}")
    best_umi = None
    best_errors = float('inf')
    best_trim_corr = float('inf')
    seq = remaining

    if len(seq) > UMI_LENGTH:
        for start in range(len(seq) - (UMI_LENGTH-1)):
            window = seq[start:start + UMI_LENGTH]
            #print(f"\nTrying window {start}-{start + 35}: {window}")
            umi_candidate, errors, trim_corr, err_base_counts, err_event_counts  = correct_and_validate(window, seq[start + UMI_LENGTH:], max_errors, max_total_errors)
            if umi_candidate:
                if errors < best_errors or (errors == best_errors and trim_corr < best_trim_corr):
                    best_umi = umi_candidate
                    best_errors = errors
                    best_trim_corr = trim_corr
                    best_base_counts = err_base_counts
                    best_event_counts = err_event_counts
    else:
        umi_candidate, errors, trim_corr, err_base_counts, err_event_counts = correct_and_validate(seq, "", max_errors, max_total_errors)
        if umi_candidate:
            best_umi = umi_candidate
            best_errors = errors
            best_trim_corr = trim_corr
            best_base_counts = err_base_counts
            best_event_counts = err_event_counts

    if best_umi:
        #print(f"\n✅ Final UMI: {best_umi} (errors: {best_errors}, trimers corrected: {best_trim_corr})")
        #return best_umi, best_errors, best_trim_corr, best_base_counts, best_event_counts
        return best_umi, best_errors, remaining, best_trim_corr, best_base_counts, best_event_counts

    else:
        #print("\n❌ No valid UMI found.")
        #return None, 0, 0, 0, 0
        return None, 0, remaining, 0, 0, 0

# ------------------ Argument Parsing ------------------

#outfile = out_fastq = iotools.open_file(args.outname, "w")
#log =  iotools.open_file(args.outname + ".log","w")
infile = args.infile
errors_permitted_in_umi = int(args.errors)

# ------------------ Key Params ------------------

UMI_TRIMERS = 12
UMI_LENGTH = UMI_TRIMERS * 3
ERRORS_ALLOWED = errors_permitted_in_umi
match_score = 2 # for parasail fuzzy match
mismatch_score = -1 # for parasail fuzzy match
matrix = parasail.matrix_create("ACGT", match_score, mismatch_score)
VALIDATION_THRESHOLD = (UMI_LENGTH * match_score) - (ERRORS_ALLOWED * (match_score - mismatch_score))

# ------------------ Adapter Pattern ------------------

adapter_pattern = "GTACTCTGCGTTGATACCACTGCTT"
error = 0.35 ## error threshold for adapter matching
primer_length = len(adapter_pattern)
min_score_threshold = (1 - error) * primer_length * match_score

# ------------------ Counters ------------------

total_reads = 0
adapter_matched = 0
umi_extracted = 0
adapter_missing = 0
umi_failed = 0
umi_candidate_lengths = Counter()
short_count = 0
long_count = 0
error_count = Counter()
G_count = 0
umi_score_counter = Counter()
low_score = 0
error_type_counter = Counter({'ins': 0, 'del': 0, 'subs': 0})


# ------------------ Main Processing ------------------
with pysam.FastxFile(infile) as fh, \
     iotools.open_file(args.outname, "w") as out_fastq:

    for n, record in enumerate(fh, 1):
        #if n > args.max_reads:
        #    break
        total_reads += 1

        seq_nano = record.sequence[-200:]
        result = parasail.sw_trace_striped_16(adapter_pattern, seq_nano, 1, 1, matrix)

        if result.score >= min_score_threshold:
            adapter_matched += 1
            aligned_length = len(result.traceback.ref.replace("-", ""))
            match_pos = result.end_ref - aligned_length + 1
            umi_start = match_pos - 60
            if result.traceback.ref[0] != "G":
                G_count += 1
                match_pos+=2

            if umi_start >= 0:
                umi_window = seq_nano[umi_start:match_pos]
                umi, total_errors, used_seq, trimers_corrected, bases_correct_error_types, trimers_corrected_error_types  = extract_umi_from_polyA_region(umi_window)
                
                if used_seq:
                    seq_len = len(used_seq)
                if seq_len < 34:
                    short_count += 1
                elif seq_len > 39:
                    long_count += 1
                else:
                    umi_candidate_lengths[seq_len] += 1

                # if len(used_seq) < 34 or len(used_seq) > 39:
                #     outlier_log.write(
                #         f"Read: {record.name}\n"
                #         f"seq_nano: {seq_nano}\n"
                #         f"used_seq: {used_seq} (length: {len(used_seq)})\n"
                #         f"{'-'*40}\n"
                #      )

                umi_pass = False
                if umi:
                    umi_validation = parasail.sw_trace_striped_16(umi, seq_nano, 1, 1, matrix)
                    umi_alignment_score = umi_validation.score
                    umi_score_counter[umi_alignment_score] += 1
                    error_count[total_errors] += 1
                    for error_type, count in bases_correct_error_types.items():
                        error_type_counter[error_type] += count  # Increment the count for the respective error type
                    if umi_alignment_score >= VALIDATION_THRESHOLD:
                      umi_pass = True
                    else:
                      low_score += 1
                if umi_pass:
                    umi_extracted += 1
                    new_name = f"{record.name}_{umi}"
                    out_fastq.write(f"@{new_name}\n{record.sequence}\n+\n{record.quality}\n")
                    #umi_out.write(f"{record.name}: {umi}: {umi_alignment_score}\n")
                    #debug_tsv.write(f"{record.name}\t{umi}\t{total_errors}\t{used_seq}\n")
                else:
                    umi_failed += 1
                    #no_umi_out.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")
                    pass
            else:
                umi_failed += 1
                #no_umi_out.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")
                pass
        else:
            adapter_missing += 1
            #no_adapter_logger.info(f"Read {n} - Alignment Score: {result.score}")
            #no_adapter_out.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")
            pass
    
with iotools.open_file(args.outname + ".log", "w") as log:
    log.write(f"The number of total reads: {total_reads}\n")
    log.write(f"The number of total reads with a polyA UMI: {umi_extracted}\n")
    log.write(f"The number of total recovered percent is: {(umi_extracted/total_reads)*100}\n")
    
    log.write("\nUMI Candidate Sequence Length Histogram\n")
    log.write("=======================================\n")
    log.write(f"Short (<30bp): {short_count}\n")
    for length, count in sorted(umi_candidate_lengths.items()):
        log.write(f"Length {length}: {count} reads\n")
    log.write(f"Long (>45bp): {long_count}\n")
    
    log.write("\nError-counts histogram (total errors):\n")
    log.write("-------------------------------------------------\n")
    for total_errors, count in sorted(error_count.items()):
        log.write(f"Total errors {total_errors}: {count} reads\n")
    
    log.write("\nUmi-alignment histogram (total score):\n")
    log.write("-------------------------------------------------\n")
    for score, count in sorted(umi_score_counter.items()):
        log.write(f"UMI score {score}: {count} reads\n")
    log.write(f"\nUmi-alignment score < {VALIDATION_THRESHOLD}: ({low_score})\n")
    
    
    log.write("\nError Type Histogram (Insertions, Deletions, Substitutions):\n")
    log.write("---------------------------------------------------------------\n")
    for error_type, count in sorted(error_type_counter.items()):
        log.write(f"{error_type.capitalize()}: {count} reads\n")
    

#log.close()
#outfile.close()
