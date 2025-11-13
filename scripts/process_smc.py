#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Processes SMC output files in a folder to generate weighted sequence logos
for a specified nucleotide/amino acid region, for a selected range of sequences within blocks.
Also identifies and outputs the single most frequent sequence (MAP estimate) and the
marginal mode (MM) sequence for each of the selected target indices.

Reads multiple files, extracts log weights from headers, normalizes weights
using LogSumExp, and generates weighted sequence logos for both nucleotide
and translated amino acid sequences within the specified region and selected block indices.
"""

import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import sys
import os
import argparse
import math
import numpy as np
from pathlib import Path
import time
import traceback


# Define standard nucleotides and amino acids
nucleotides = ["A", "C", "G", "T"]
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Define Log Negative Infinity
LOG_NEG_INFINITY = -float('inf')

# --- Helper Functions for LogSumExp ---

def log_add(log_a, log_b):
    """ Calculates log(exp(a) + exp(b)) in a numerically stable way. """
    if log_a == LOG_NEG_INFINITY and log_b == LOG_NEG_INFINITY: return LOG_NEG_INFINITY
    if log_a < log_b: log_a, log_b = log_b, log_a
    if log_b == LOG_NEG_INFINITY: return log_a
    try: return log_a + math.log1p(math.exp(log_b - log_a))
    except OverflowError: print("Warning: Potential overflow in log_add."); return log_a

def log_sum_exp(log_values):
    """ Calculates log(sum(exp(v) for v in log_values)) stably using pairwise log_add. """
    valid_log_values = [lv for lv in log_values if lv is not None and lv != LOG_NEG_INFINITY]
    if not valid_log_values: return LOG_NEG_INFINITY
    current_log_sum = valid_log_values[0]
    for i in range(1, len(valid_log_values)): current_log_sum = log_add(current_log_sum, valid_log_values[i])
    return current_log_sum

# --- Argument Parsing ---

def parse_arguments():
    """ Parses command-line arguments. """
    parser = argparse.ArgumentParser(description="Process SMC output files for weighted sequence logos and identify highest frequency sequences.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder_name", type=str, help="Folder containing .out files.")
    parser.add_argument("start_index", type=int, help="Start index for files (e.g., 1).")
    parser.add_argument("end_index", type=int, help="End index for files (e.g., 45).")
    parser.add_argument("--aa_start", type=int, default=99, help="Start AA position (1-based).")
    parser.add_argument("--aa_end", type=int, default=111, help="End AA position (1-based).")
    parser.add_argument("-o", "--output_prefix", type=str, default="logo", help="Output file prefix.")
    # New argument for specific target indices
    parser.add_argument("--target_indices", type=int, nargs='+', default=None,
                        help="Space-separated list of 1-based indices of sequences in each block for logo generation (e.g., 1 5 10). Overrides --logo_start_index and --logo_end_index.")
    # Old arguments for range, will be used if --target_indices is not set
    parser.add_argument("--logo_start_index", type=int, default=1, help="1-based start index of sequence in block for logo generation (used if --target_indices is not provided).")
    parser.add_argument("--logo_end_index", type=int, default=None, help="1-based end index of sequence in block for logo generation (inclusive, used if --target_indices is not provided). Defaults to all sequences from logo_start_index.")

    args = parser.parse_args()
    if args.start_index > args.end_index: parser.error("start_index > end_index.")
    if args.aa_start > args.aa_end: parser.error("aa_start > aa_end.")
    if args.aa_start < 1: parser.error("aa_start must be >= 1.")
    
    if args.target_indices:
        if not all(idx >= 1 for idx in args.target_indices):
            parser.error("--target_indices must all be >= 1.")
    else: # Validate range arguments only if target_indices is not used
        if args.logo_start_index < 1: parser.error("logo_start_index must be >= 1.")
        if args.logo_end_index is not None and args.logo_end_index < args.logo_start_index:
            parser.error("logo_end_index must be >= logo_start_index.")
    return args

# --- File Parsing Functions ---
def parse_header_and_get_log_weight(file_handle):
    """ Reads header, returns log_weight and sequence start position. """
    header_lines = []; last_line_read = None; log_weight = None; seq_start_pos = file_handle.tell()
    while True:
        current_pos = file_handle.tell(); line = file_handle.readline()
        if not line: break
        stripped = line.strip(); is_seq = False
        if stripped:
            pot_seq_chars = set(nucleotides + ['N', 'U'])
            if all(c.upper() in pot_seq_chars for c in stripped):
                if len(stripped) > 10: is_seq = True
            elif '(' in stripped or 'ESS' in stripped: is_seq = False # Heuristic for some header lines
            else:
                try:
                    float(stripped) # Check if it's a number
                    if all(c in '0123456789.-+Ee' for c in stripped): is_seq = False # If purely numerical, not sequence
                    else: is_seq = True # Numerical but with other chars, likely sequence
                except ValueError:
                    if stripped: is_seq = True # Not a number, non-empty, assume sequence or other header
        if is_seq: file_handle.seek(current_pos); seq_start_pos = current_pos; break
        else:
            if stripped: header_lines.append(stripped); last_line_read = stripped
            seq_start_pos = file_handle.tell() # Update in case this is the end of header/file
    if last_line_read:
        try: log_weight = float(last_line_read.split()[0]); return log_weight, seq_start_pos
        except (ValueError, IndexError): print(f"Warn: Bad log weight line: '{last_line_read}'."); return None, seq_start_pos
    else: print("Warn: No header lines found."); return None, 0


def validate_block(block, expected_len, line_num):
    """ Validates block content. """
    if not block: return False, 0, 0
    first_len = len(block[0]); block_len = len(block)
    if first_len != expected_len: return False, block_len, first_len
    if not all(len(seq) == first_len for seq in block): return False, block_len, first_len
    return True, block_len, first_len

def parse_sequence_blocks(file_handle, seq_start_pos, start_nt, end_nt, expected_block_size=None):
    """ Reads sequence blocks starting from seq_start_pos. """
    try: file_handle.seek(seq_start_pos)
    except ValueError: print(f"Error: Invalid seek pos {seq_start_pos}."); return [], None
    max_len = end_nt; extract_len = end_nt - start_nt; nt_set = set(nucleotides)
    blocks = []; reads = 0; skip_short = 0; skip_incon = 0; file_block_size = None
    current_block = []; line_num = 0
    while True:
        try: line = file_handle.readline()
        except Exception as e: print(f"Error reading line: {e}"); break
        if not line: break
        line_num += 1; line = line.strip()
        if line:
            reads += 1
            if len(line) < max_len: skip_short += 1; continue
            sub_seq = line[start_nt:end_nt]
            if len(sub_seq) != extract_len: skip_incon += 1; continue
            if not all(c.upper() in nt_set for c in sub_seq): skip_incon += 1; continue
            current_block.append(sub_seq.upper()) # Store uppercase
        else:
            if current_block:
                ok, blen, slen = validate_block(current_block, extract_len, line_num)
                if ok:
                    if file_block_size is None: file_block_size = blen
                    if expected_block_size is not None and file_block_size != expected_block_size: skip_incon += len(current_block)
                    elif blen == file_block_size: blocks.append(current_block)
                    else: skip_incon += len(current_block)
                else: skip_incon += len(current_block)
                current_block = []
    if current_block: # Process last block if file doesn't end with blank line
        ok, blen, slen = validate_block(current_block, extract_len, line_num)
        if ok:
            if file_block_size is None: file_block_size = blen
            if expected_block_size is not None and file_block_size != expected_block_size: skip_incon += len(current_block)
            elif blen == file_block_size: blocks.append(current_block)
            else: skip_incon += len(current_block)
        else: skip_incon += len(current_block)
    # if skip_short > 0 or skip_incon > 0:
    #     print(f"  parse_sequence_blocks: Total lines read={reads}, skipped short={skip_short}, inconsistent/invalid char={skip_incon}")
    return blocks, file_block_size

# --- PFM Handling Functions ---

def initialize_pfms(n_logos_to_generate, symbols, seq_len):
    """ Creates a list of empty PFM DataFrames. """
    if n_logos_to_generate is None or n_logos_to_generate <= 0:
        print("Info: No logos to generate based on selected indices or block size.")
        return []
    if seq_len <= 0: print(f"Error: Invalid seq length ({seq_len}). Cannot initialize PFMs."); return None # Critical error
    print(f"Initializing {n_logos_to_generate} PFMs for {len(symbols)} symbols (length {seq_len}).")
    return [pd.DataFrame(0.0, index=symbols, columns=range(seq_len)) for _ in range(n_logos_to_generate)]


def update_pfms(pfm_list_df, list_of_blocks_from_file, weight, valid_symbols, seq_length,
                is_aa_data=False, nt_indices_to_extract_0based=None, expected_raw_block_size_for_nt=None):
    """
    Updates PFM counts.
    If is_aa_data=True, list_of_blocks_from_file is a list of lists, where each inner list
                     contains AA sequences (or None) for each PFM.
    If is_aa_data=False (NT), list_of_blocks_from_file is a list of raw full NT blocks,
                               and nt_indices_to_extract_0based specifies which sequences to use.
    """
    if not pfm_list_df: return 0
    if not list_of_blocks_from_file: return 0

    num_pfms_to_update = len(pfm_list_df)
    if num_pfms_to_update == 0: return 0

    pfm_shape_cols = pfm_list_df[0].shape[1]
    if seq_length != pfm_shape_cols:
        print(f"Internal Error: PFM sequence length mismatch ({seq_length} vs {pfm_shape_cols}). Cannot update.")
        return 0

    num_symbols = len(valid_symbols)
    symbol_to_int = {symbol: i for i, symbol in enumerate(valid_symbols)}
    temp_counts_np = [np.zeros((num_symbols, seq_length), dtype=np.float64) for _ in range(num_pfms_to_update)]
    sequences_processed_for_pfm_update = 0

    for block_data in list_of_blocks_from_file:
        if is_aa_data: # AA data: block_data is a list of AA sequences or None (one such list per "raw block")
            if len(block_data) != num_pfms_to_update:
                print(f"  Warn: AA data block has {len(block_data)} items, expected {num_pfms_to_update} for PFMs. Skipping this AA data block.")
                continue
            for pfm_idx, seq in enumerate(block_data):
                if seq is None:  # Skip failed translations or inconsistent lengths
                    continue
                # Length check is crucial here, should match the PFM's seq_length
                if len(seq) != seq_length:
                    print(f"  Warn: AA Seq len mismatch {len(seq)} vs {seq_length} for PFM {pfm_idx}. Skipping seq.")
                    continue
                for pos, symbol in enumerate(seq):
                    upper_symbol = symbol.upper()
                    if upper_symbol in symbol_to_int:
                        row_idx = symbol_to_int[upper_symbol]
                        temp_counts_np[pfm_idx][row_idx, pos] += weight
                sequences_processed_for_pfm_update += 1
        else: # NT data: block_data is a raw full NT block (list of NT sequences)
            if expected_raw_block_size_for_nt is not None and len(block_data) != expected_raw_block_size_for_nt:
                print(f"  Warn: NT raw block has {len(block_data)} seqs, expected {expected_raw_block_size_for_nt}. Skipping block.")
                continue
            if nt_indices_to_extract_0based is None:
                print("  Error: nt_indices_to_extract_0based is None for NT data processing. Skipping block.")
                continue

            for pfm_idx, raw_seq_idx in enumerate(nt_indices_to_extract_0based):
                if raw_seq_idx >= len(block_data):
                    print(f"  Warn: NT raw sequence index {raw_seq_idx} out of bounds for block of size {len(block_data)}. Skipping for PFM {pfm_idx}.")
                    continue
                seq = block_data[raw_seq_idx] # This is the pre-sliced NT sequence
                if len(seq) != seq_length:
                    print(f"  Warn: NT Seq len mismatch {len(seq)} vs {seq_length} for PFM {pfm_idx} (raw index {raw_seq_idx}). Skipping seq.")
                    continue
                for pos, symbol in enumerate(seq):
                    upper_symbol = symbol.upper() # Should already be upper from parse_sequence_blocks
                    if upper_symbol in symbol_to_int: # nucleotides
                        row_idx = symbol_to_int[upper_symbol]
                        temp_counts_np[pfm_idx][row_idx, pos] += weight
                sequences_processed_for_pfm_update +=1
    try:
        for pfm_idx in range(num_pfms_to_update):
            pfm_list_df[pfm_idx] += temp_counts_np[pfm_idx]
    except Exception as e:
        print(f"Error adding numpy counts to DataFrame PFM {pfm_idx}: {e}")
        traceback.print_exc()
    return sequences_processed_for_pfm_update


def generate_probability_matrices(pfm_list):
    """ Converts weighted PFMs to probability matrices. """
    if not pfm_list: print("Warn: No PFMs to generate probabilities from."); return []
    prob_mats = []; print(f"Calculating probabilities from {len(pfm_list)} PFMs...")
    for i, pfm in enumerate(pfm_list):
        if pfm.empty: prob_mats.append(pd.DataFrame()); continue
        tot_w = pfm.sum(axis=0); prob_m = pfm.copy()
        for pos in prob_m.columns:
            tw = tot_w[pos]
            if tw > 1e-9: prob_m[pos] = prob_m[pos] / tw
            else: prob_m[pos] = 0.0 # Avoid division by zero; all symbols at this position have zero weight
        prob_mats.append(prob_m.T) # Transpose for logomaker (symbols as columns, positions as rows)
    print(f"Generated {len(prob_mats)} probability matrices.")
    return prob_mats

# --- Translation Function ---
def translate_sequence_blocks(nt_blocks, selected_indices_0based_in_block, expected_aa_len_from_args, expected_raw_block_size):
    """
    Translates specific nucleotide sequences (selected by selected_indices_0based_in_block)
    from nt_blocks to amino acid sequences.
    Returns a list of "AA blocks", where each AA block is a list of translated sequences or None.
    Also returns the determined actual_aa_len for consistent processing.
    """
    translated_aa_blocks = []
    total_translation_attempts = 0
    translation_errors = 0

    if not nt_blocks: return [], expected_aa_len_from_args # Return arg-based len if no NT blocks

    # Determine actual_aa_len from the first successfully translated selected sequence that matches expected_aa_len_from_args if possible
    actual_aa_len_for_pfm = None
    first_valid_aa_len_found = False

    # Try to establish a definitive actual_aa_len using expected_aa_len_from_args as a guide
    if expected_aa_len_from_args is not None and expected_aa_len_from_args > 0:
        for nt_block_for_len_check in nt_blocks:
            if len(nt_block_for_len_check) != expected_raw_block_size: continue
            for idx_in_block in selected_indices_0based_in_block:
                if idx_in_block < len(nt_block_for_len_check):
                    nt_seq_for_len = nt_block_for_len_check[idx_in_block]
                    if len(nt_seq_for_len) % 3 == 0 and len(nt_seq_for_len) > 0:
                        try:
                            temp_trans = str(Seq(nt_seq_for_len).translate(to_stop=True)) # Use to_stop=True for length check
                            if len(temp_trans) == expected_aa_len_from_args:
                                actual_aa_len_for_pfm = expected_aa_len_from_args
                                print(f"  Confirmed actual AA length for PFM: {actual_aa_len_for_pfm} (matches expected from args) from selected NT seq (index {idx_in_block}).")
                                first_valid_aa_len_found = True
                                break
                        except Exception: pass
            if first_valid_aa_len_found: break
    
    # If not confirmed via expected_aa_len_from_args, try to determine from any translation
    if not first_valid_aa_len_found:
        for nt_block_for_len_check in nt_blocks:
            if len(nt_block_for_len_check) != expected_raw_block_size: continue
            for idx_in_block in selected_indices_0based_in_block:
                if idx_in_block < len(nt_block_for_len_check):
                    nt_seq_for_len = nt_block_for_len_check[idx_in_block]
                    if len(nt_seq_for_len) % 3 == 0 and len(nt_seq_for_len) > 0:
                        try:
                            temp_trans_len = len(str(Seq(nt_seq_for_len).translate(to_stop=True)))
                            if temp_trans_len > 0:
                                actual_aa_len_for_pfm = temp_trans_len
                                print(f"  Determined actual AA length for PFM: {actual_aa_len_for_pfm} from selected NT seq (index {idx_in_block}). May differ from arg-based length.")
                                first_valid_aa_len_found = True
                                break
                        except Exception: pass
            if first_valid_aa_len_found: break

    if not first_valid_aa_len_found:
        if expected_aa_len_from_args is not None and expected_aa_len_from_args > 0:
            actual_aa_len_for_pfm = expected_aa_len_from_args # Fallback to arg-based if primary checks failed
            print(f"  Warning: Could not determine AA length from translated sequences. Using expected_aa_len_from_args: {actual_aa_len_for_pfm}")
        else:
            print(f"  Critical Warning: Could not determine a valid AA length for translation. AA processing will likely fail or be inconsistent.")
            # actual_aa_len_for_pfm remains None, subsequent steps should handle this.

    for blk_idx, nt_raw_block in enumerate(nt_blocks):
        if len(nt_raw_block) != expected_raw_block_size:
            print(f"  Warn: NT Raw Block {blk_idx} has {len(nt_raw_block)} seqs, expected {expected_raw_block_size}. Skipping this raw block for translation.")
            # Need to decide if we append a list of Nones or skip this block entirely from translated_aa_blocks.
            # For consistency with PFM updates, if a raw block is skipped, its corresponding AA "block" shouldn't be created.
            continue

        current_translated_aa_selection = [] # Holds AA seqs for one raw NT block, for the selected indices
        for target_seq_0based_idx in selected_indices_0based_in_block:
            total_translation_attempts += 1
            aa_seq_to_add = None # Default to None

            if target_seq_0based_idx >= len(nt_raw_block):
                # This case should ideally not happen if nt_raw_block conforms to expected_raw_block_size
                # and selected_indices_0based_in_block are within this size.
                # However, good to have a check.
                print(f"  Warn: Target NT index {target_seq_0based_idx} out of bounds for raw block {blk_idx} of size {len(nt_raw_block)}. Adding None for this AA sequence.")
                translation_errors += 1
            else:
                nt_seq_to_translate = nt_raw_block[target_seq_0based_idx]
                if len(nt_seq_to_translate) % 3 != 0:
                    print(f"  Warn: NT Seq (len {len(nt_seq_to_translate)}) for AA target index {selected_indices_0based_in_block.index(target_seq_0based_idx)+1} (raw NT idx {target_seq_0based_idx}) in raw block {blk_idx} not div by 3. Adding None.")
                    translation_errors += 1
                elif not nt_seq_to_translate: # Empty sequence
                    print(f"  Warn: NT Seq for AA target index {selected_indices_0based_in_block.index(target_seq_0based_idx)+1} (raw NT idx {target_seq_0based_idx}) in raw block {blk_idx} is empty. Adding None.")
                    translation_errors += 1
                else:
                    try:
                        translated_aa_seq = str(Seq(nt_seq_to_translate).translate(to_stop=True)) # Translate with stop
                        if actual_aa_len_for_pfm is not None and len(translated_aa_seq) != actual_aa_len_for_pfm:
                            # This condition means the current translation's length does not match the *established* AA length for PFMs
                            print(f"  Warn: Translated AA seq (len {len(translated_aa_seq)}) from NT (raw idx {target_seq_0based_idx}) for AA target {selected_indices_0based_in_block.index(target_seq_0based_idx)+1} in raw block {blk_idx} differs from established PFM AA length ({actual_aa_len_for_pfm}). Adding None.")
                            translation_errors += 1
                        elif actual_aa_len_for_pfm is None and len(translated_aa_seq) > 0: # This is the first valid one, establish length
                            # This case should ideally be caught by the initial length determination,
                            # but as a fallback if the first file had no valid translations.
                            actual_aa_len_for_pfm = len(translated_aa_seq)
                            print(f"  Established PFM AA length on the fly: {actual_aa_len_for_pfm} from NT (raw idx {target_seq_0based_idx}), raw block {blk_idx}")
                            aa_seq_to_add = translated_aa_seq.upper()
                        elif actual_aa_len_for_pfm is None and len(translated_aa_seq) == 0 : # Could not establish length
                             print(f"  Warn: Translated AA seq from NT (raw idx {target_seq_0based_idx}) for AA target {selected_indices_0based_in_block.index(target_seq_0based_idx)+1} in raw block {blk_idx} is empty, and no PFM AA length established. Adding None.")
                             translation_errors +=1
                        else: # Length matches or was just established (and not None)
                            aa_seq_to_add = translated_aa_seq.upper()
                    except Exception as e:
                        print(f"  Error translating NT Seq from NT (raw idx {target_seq_0based_idx}) for AA target {selected_indices_0based_in_block.index(target_seq_0based_idx)+1} in raw block {blk_idx}: {e}. Adding None.")
                        translation_errors += 1
            current_translated_aa_selection.append(aa_seq_to_add)
        
        translated_aa_blocks.append(current_translated_aa_selection)

    if translation_errors > 0:
        print(f"Translation: {translation_errors}/{total_translation_attempts} selected sequences had errors/issues or inconsistent lengths.")
    
    return translated_aa_blocks, actual_aa_len_for_pfm


# --- Plotting Function ---
def plot_sequence_logos(prob_mats, base_filename,
                        plot_title_indices_1based, # List of 1-based indices for subplot titles
                        is_aa=False, start_lbl_for_aa_axis=None,
                        max_plots_per_figure=10):
    """ Generates and saves sequence logo plots in PDF format with specific styling, chunking large sets. """
    if not prob_mats:
        # base_title derived from base_filename for logging
        base_title = Path(base_filename).stem
        print(f"No data to plot for {base_title}.")
        return

    num_total_mats = len(prob_mats)
    if num_total_mats == 0 :
        base_title = Path(base_filename).stem
        print(f"No probability matrices to plot for {base_title}.")
        return
    
    if len(plot_title_indices_1based) != num_total_mats:
        print(f"  Error: Mismatch between number of probability matrices ({num_total_mats}) and title indices ({len(plot_title_indices_1based)}). Cannot generate plot titles correctly.")
        # Fallback: generate generic titles or skip. For now, let it proceed, titles might be off.
        # Or, more strictly: return

    for chunk_num, chunk_start_idx in enumerate(range(0, num_total_mats, max_plots_per_figure)):
        chunk_end_idx = min(chunk_start_idx + max_plots_per_figure, num_total_mats)
        current_chunk_prob_mats = prob_mats[chunk_start_idx:chunk_end_idx]
        n_plots_in_chunk = len(current_chunk_prob_mats)

        if n_plots_in_chunk == 0: continue

        path_obj = Path(base_filename)
        if num_total_mats > max_plots_per_figure:
            chunk_filename_str = f"{path_obj.stem}_part{chunk_num + 1}{path_obj.suffix}"
        else:
            chunk_filename_str = base_filename
        
        chunk_file_path = path_obj.parent / chunk_filename_str
        
        seq_len_for_fig_calc = 20 # Default
        first_valid_mat_in_chunk = next((m for m in current_chunk_prob_mats if not m.empty), None)
        if first_valid_mat_in_chunk is not None:
            seq_len_for_fig_calc = first_valid_mat_in_chunk.shape[0] # pfm is transposed, so shape[0] is length

        fig_w = max(12, seq_len_for_fig_calc * 0.7)
        fig_h = max(4, n_plots_in_chunk * 3.5) # Increased height per plot slightly

        fig, axes = plt.subplots(nrows=n_plots_in_chunk, ncols=1, figsize=(fig_w, fig_h), squeeze=False)

        for i_in_chunk, pfm_df_transposed in enumerate(current_chunk_prob_mats):
            ax = axes[i_in_chunk, 0]
            original_pfm_idx = chunk_start_idx + i_in_chunk
            
            current_plot_actual_sequence_index = plot_title_indices_1based[original_pfm_idx]
            plot_title_text = f"Marginal Posterior - Target Index {current_plot_actual_sequence_index}"


            if pfm_df_transposed.empty:
                ax.set_title(f"{plot_title_text}: Empty", fontsize=16, weight='bold') # Adjusted font size
                ax.text(0.5,0.5,'No data',ha='center',va='center')
                ax.set_xticks([]); ax.set_yticks([])
                for spine in ax.spines.values(): spine.set_visible(False)
                continue

            pfp = pfm_df_transposed # Already transposed (positions as rows, symbols as columns)
            if not isinstance(pfp.index, pd.RangeIndex): # Ensure index is 0, 1, 2... for logomaker
                 pfp = pfp.reset_index(drop=True)

            try:
                logo_color_scheme = 'hydrophobicity' if is_aa else 'classic'
                logo = logomaker.Logo(pfp,
                                      ax=ax,
                                      font_name='Arial Rounded MT Bold',
                                      color_scheme=logo_color_scheme)
                
                ax.set_title(plot_title_text, fontsize=16, weight='bold') # Adjusted font size
                ax.set_ylabel("Posterior Probability", fontsize=14) # Adjusted font size

                logo.style_spines(visible=False)
                logo.style_spines(spines=['left', 'bottom'], visible=True, linewidth=1.5)

                ax.tick_params(axis='both', which='major', labelsize=10) # Adjusted font size

                max_stack_height = 1.0 # For probability matrices
                # if not pfp.empty: # This check is not strictly needed for probability matrices as sum is ~1
                #     pass

                ylim_max = max(max_stack_height * 1.05, 0.1)
                ax.set_ylim([0, ylim_max])

                logo.style_xticks(rotation=0, fmt='%d', anchor=0)

                current_tick_values = ax.get_xticks()
                visible_ticks_float = [tick for tick in current_tick_values if tick >= 0 and tick < len(pfp)]
                visible_ticks = sorted(list(set(map(lambda x: int(round(x)), visible_ticks_float))))
                visible_ticks = [t for t in visible_ticks if t >=0 and t < len(pfp)]


                if is_aa and start_lbl_for_aa_axis is not None:
                    new_labels = []
                    actual_ticks_to_label = []
                    for tick_val in visible_ticks:
                        actual_aa_pos = start_lbl_for_aa_axis + tick_val
                        new_labels.append(str(actual_aa_pos))
                        actual_ticks_to_label.append(tick_val)
                    ax.set_xticks(actual_ticks_to_label)
                    ax.set_xticklabels(new_labels)
                    ax.set_xlabel(f"Amino Acid Position (Original: {start_lbl_for_aa_axis}+)", fontsize=14)
                else: # NT
                    new_labels = [str(tick_val + 1) for tick_val in visible_ticks]
                    ax.set_xticks(visible_ticks)
                    ax.set_xticklabels(new_labels)
                    ax.set_xlabel("Nucleotide Position (in Slice)", fontsize=14)
                
                num_total_positions_in_logo = len(pfp)
                rotation_angle = 0
                if num_total_positions_in_logo > 25 and len(visible_ticks) > 15:
                     rotation_angle = 45
                
                if rotation_angle > 0:
                    plt.setp(ax.get_xticklabels(), rotation=rotation_angle, ha="right", rotation_mode="anchor")
                else:
                    plt.setp(ax.get_xticklabels(), ha="center")

            except Exception as e:
                print(f"\nError during logo plotting for {plot_title_text}: {e}")
                traceback.print_exc()
                ax.set_title(f"{plot_title_text}: Plot Error", fontsize=16, weight='bold')
                ax.text(0.5,0.5,'Plot Error',ha='center',va='center',color='red')

        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        try:
            chunk_file_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(chunk_file_path, dpi=300, bbox_inches='tight')
            print(f"\nLogo chunk saved: '{chunk_file_path}'.")
        except Exception as e:
            print(f"Error saving '{chunk_file_path}': {e}")
        plt.close(fig)

# --- Function to find and print highest frequency sequences ---
def find_and_print_highest_frequency_sequences(sequence_overall_weights_list, sequence_type_str, plot_title_indices_1based, aa_start_arg=None):
    """
    Finds and prints the sequence with the highest accumulated normalized weight for each target index.
    This corresponds to the Maximum a Posteriori (MAP) estimate.
    """
    print(f"\n--- Highest Frequency (MAP) {sequence_type_str} Sequences (Overall Accumulated Weight) ---")
    if not sequence_overall_weights_list:
        print(f"No {sequence_type_str} sequence frequency data to report.")
        return

    found_any = False
    for pfm_idx, seq_weights_dict in enumerate(sequence_overall_weights_list):
        target_idx_1_based = plot_title_indices_1based[pfm_idx]
        
        title = f"For Target Index {target_idx_1_based}"
        if sequence_type_str == "Amino Acid" and aa_start_arg is not None:
            title += f" (Input AA Start Ref: {aa_start_arg})"

        if not seq_weights_dict:
            print(f"\n{title}:\n  No sequences recorded for this target index.")
            continue

        # Find the sequence with the maximum weight
        # Sort by weight (desc) then by sequence (asc) for tie-breaking
        sorted_sequences = sorted(seq_weights_dict.items(), key=lambda item: (-item[1], item[0]))
        
        if not sorted_sequences: # Should not happen if seq_weights_dict is not empty
            print(f"\n{title}:\n  No weighted sequences found after sorting.")
            continue

        highest_freq_seq, highest_weight = sorted_sequences[0]
        
        print(f"\n{title}:")
        print(f"  Sequence: {highest_freq_seq}")
        print(f"  Accumulated Normalized Weight: {highest_weight:.6e}")
        if len(sorted_sequences) > 1:
            second_freq_seq, second_weight = sorted_sequences[1]
            # Check for ties or very close weights (e.g., within machine epsilon for typical float operations)
            if abs(highest_weight - second_weight) < max(1e-9 * max(abs(highest_weight), abs(second_weight)), 1e-12) :
                print(f"  Note: Other sequences have similar high weights (e.g., '{second_freq_seq}' with weight {second_weight:.6e}).")
        found_any = True

    if not found_any: # This means all seq_weights_dict were empty
        print(f"No specific highest frequency {sequence_type_str} sequences found for any target index.")


# --- Function to find and print marginal mode sequences ---
def find_and_print_marginal_mode_sequences(pfm_list, sequence_type_str, plot_title_indices_1based, aa_start_arg=None):
    """
    Finds and prints the Marginal Mode (MM) sequence for each target index.
    The MM sequence is constructed by taking the most probable character at each position independently.
    """
    print(f"\n--- Marginal Mode (MM) {sequence_type_str} Sequences ---")
    if not pfm_list or not any(not pfm.empty for pfm in pfm_list):
        print(f"No {sequence_type_str} PFM data to report for MM estimate.")
        return

    for pfm_idx, pfm_df in enumerate(pfm_list):
        target_idx_1_based = plot_title_indices_1based[pfm_idx]

        title = f"For Target Index {target_idx_1_based}"
        if sequence_type_str == "Amino Acid" and aa_start_arg is not None:
            title += f" (Input AA Start Ref: {aa_start_arg})"
        
        print(f"\n{title}:")

        if pfm_df.empty:
            print("  No PFM data recorded for this target index.")
            continue

        try:
            # For each column (position), find the index (character) with the maximum value.
            # This gives a pandas Series of the most probable characters at each position.
            mm_sequence_chars = pfm_df.idxmax(axis=0)
            
            # Join the characters to form the final MM sequence string.
            mm_sequence = "".join(mm_sequence_chars)
            
            print(f"  Sequence: {mm_sequence}")

        except Exception as e:
            print(f"  Could not calculate MM sequence due to an error: {e}")


# --- Main Execution Logic ---
def main():
    args = parse_arguments()
    folder_path = Path(args.folder_name)
    if not folder_path.is_dir(): print(f"Error: Folder '{args.folder_name}' not found."); sys.exit(1)

    start_aa_arg, end_aa_arg = args.aa_start, args.aa_end
    start_nt_slice_0based, end_nt_slice_0based_exclusive = (start_aa_arg - 1) * 3, end_aa_arg * 3
    nt_region_len = end_nt_slice_0based_exclusive - start_nt_slice_0based
    
    aa_region_len_from_args = (end_aa_arg - start_aa_arg) + 1 if end_aa_arg >= start_aa_arg else 0

    if nt_region_len <= 0: print("Error: Invalid AA range -> non-positive nucleotide region length."); sys.exit(1)

    print(f"Processing: '{folder_path}', File Index {args.start_index}-{args.end_index}, Input AA Region: {start_aa_arg}-{end_aa_arg} (NT Slice {start_nt_slice_0based}-{end_nt_slice_0based_exclusive-1})")
    print(f"Nucleotide region length for PFM/Logo: {nt_region_len}")
    if aa_region_len_from_args > 0:
        print(f"Expected Amino Acid region length (from args): {aa_region_len_from_args}")
    else:
        print("Amino Acid processing will be skipped (aa_start/aa_end implies zero or negative length).")

    print("\n--- Pass 1: Reading headers and determining block structure ---")
    log_w_data = []; n_total_seqs_per_block = None; files_chk = 0; files_ok_hdr = 0
    for i in range(args.start_index, args.end_index + 1):
        fpath = folder_path / f"{i}.out"; files_chk += 1
        if not fpath.is_file(): print(f"Warn: File '{fpath.name}' not found, skipping."); continue
        try:
            with open(fpath, "r") as f:
                log_w, seq_start = parse_header_and_get_log_weight(f)
                if log_w is None or seq_start < 0: print(f"Warn: Invalid header/no sequences in '{fpath.name}'."); continue

                if n_total_seqs_per_block is None:
                    _, tmp_bsize = parse_sequence_blocks(f, seq_start, start_nt_slice_0based, end_nt_slice_0based_exclusive, None)
                    if tmp_bsize is not None and tmp_bsize > 0:
                        n_total_seqs_per_block = tmp_bsize
                        print(f"Determined block size (sequences per block): {n_total_seqs_per_block} (from '{fpath.name}')")
                    else:
                        print(f"Warn: Could not determine block size from '{fpath.name}'. Skipping this file for block size determination."); continue
                log_w_data.append({"path": fpath, "log_w": log_w, "start_pos": seq_start}); files_ok_hdr += 1
        except Exception as e: print(f"Error Pass 1 ({fpath.name}): {e}")

    print(f"Pass 1: Checked={files_chk}, Files with Valid Headers={files_ok_hdr}.")
    if not log_w_data : print("\nError: No files with valid headers found in Pass 1."); sys.exit(1)
    if n_total_seqs_per_block is None: print("\nError: Could not determine sequence block size."); sys.exit(1)

    selected_indices_0based_in_block = []
    plot_title_indices_1based = []

    if args.target_indices:
        valid_target_indices_1based = []
        for ti in sorted(list(set(args.target_indices))):
            if ti > n_total_seqs_per_block:
                print(f"Warning: Target index {ti} > sequences per block ({n_total_seqs_per_block}). Skipping.")
            else:
                valid_target_indices_1based.append(ti)
        
        if not valid_target_indices_1based:
            print("Error: No valid target indices. Exiting."); sys.exit(1)
        
        selected_indices_0based_in_block = [idx - 1 for idx in valid_target_indices_1based]
        plot_title_indices_1based = valid_target_indices_1based
        print(f"Will process specific sequences at 1-based indices: {plot_title_indices_1based} within each block.")
    else:
        logo_idx_start_1based_range = args.logo_start_index
        logo_idx_end_1based_range = args.logo_end_index if args.logo_end_index is not None else n_total_seqs_per_block

        if logo_idx_start_1based_range > n_total_seqs_per_block:
            print(f"Error: --logo_start_index ({logo_idx_start_1based_range}) > sequences per block ({n_total_seqs_per_block})."); sys.exit(1)
        if logo_idx_end_1based_range > n_total_seqs_per_block:
            print(f"Warning: --logo_end_index ({logo_idx_end_1based_range}) > sequences per block ({n_total_seqs_per_block}). Adjusting to {n_total_seqs_per_block}.")
            logo_idx_end_1based_range = n_total_seqs_per_block
        if logo_idx_start_1based_range > logo_idx_end_1based_range :
            print(f"Error: logo_start_index ({logo_idx_start_1based_range}) > logo_end_index ({logo_idx_end_1based_range})."); sys.exit(1)
        
        selected_indices_0based_in_block = list(range(logo_idx_start_1based_range - 1, logo_idx_end_1based_range))
        plot_title_indices_1based = list(range(logo_idx_start_1based_range, logo_idx_end_1based_range + 1))
        print(f"Will process sequences from index {logo_idx_start_1based_range} to {logo_idx_end_1based_range} (1-based) within each block.")

    num_logos_to_process = len(selected_indices_0based_in_block)
    if num_logos_to_process == 0:
        print("No sequences selected for processing. Exiting."); sys.exit(0)

    log_ws = [d["log_w"] for d in log_w_data]; log_sum_w = log_sum_exp(log_ws)
    if log_sum_w == LOG_NEG_INFINITY: print("\nError: LogSumExp of weights is -infinity. Cannot normalize."); sys.exit(1)
    print(f"\nLogSumExp of file weights: {log_sum_w:.6f}")

    overall_nt_pfms = initialize_pfms(num_logos_to_process, nucleotides, nt_region_len)
    overall_aa_pfms = []
    
    if overall_nt_pfms is None:
         print("Error initializing NT PFMs. Exiting."); sys.exit(1)
    if nt_region_len > 0 and not overall_nt_pfms and num_logos_to_process > 0 :
         print("Error: NT PFMs list empty despite valid parameters. Exiting."); sys.exit(1)

    nt_sequence_overall_weights = [{} for _ in range(num_logos_to_process)]
    aa_sequence_overall_weights = [{} for _ in range(num_logos_to_process)]

    print("\n--- Pass 2: Processing sequences, updating PFMs, and tracking sequence frequencies ---")
    files_proc_p2 = 0; total_norm_w_sum = 0.0
    actual_aa_len_for_pfm_established = None
    aa_pfms_initialized_with_correct_length = False

    for file_idx, file_data in enumerate(log_w_data):
        fpath = file_data["path"]; log_w = file_data["log_w"]; seq_start = file_data["start_pos"]
        norm_w = math.exp(log_w - log_sum_w) if log_sum_w != LOG_NEG_INFINITY else 0.0
        total_norm_w_sum += norm_w
        print(f"\n[{file_idx+1}/{len(log_w_data)}] Processing: '{fpath.name}' (Norm W: {norm_w:.3e})")

        try:
            with open(fpath, "r") as f:
                parse_start_time = time.time()
                current_raw_nt_blocks, blk_size_from_file = parse_sequence_blocks(f, seq_start, start_nt_slice_0based, end_nt_slice_0based_exclusive, n_total_seqs_per_block)
                # print(f"  parse_sequence_blocks: {parse_end_time - parse_start_time:.2f}s -> {len(current_raw_nt_blocks)} NT blocks, file_bsize {blk_size_from_file}")

                if not current_raw_nt_blocks:
                    print(f"  Skipping '{fpath.name}': No valid NT blocks parsed."); continue
                if blk_size_from_file != n_total_seqs_per_block:
                    print(f"  Skipping '{fpath.name}': Parsed block size ({blk_size_from_file}) != expected ({n_total_seqs_per_block})."); continue
                files_proc_p2 += 1

                if overall_nt_pfms and nt_region_len > 0:
                    update_pfms(overall_nt_pfms, current_raw_nt_blocks, norm_w, nucleotides, nt_region_len,
                                is_aa_data=False, nt_indices_to_extract_0based=selected_indices_0based_in_block,
                                expected_raw_block_size_for_nt=n_total_seqs_per_block)
                    for pfm_idx, target_raw_seq_idx in enumerate(selected_indices_0based_in_block):
                        for nt_block in current_raw_nt_blocks:
                            if target_raw_seq_idx < len(nt_block):
                                nt_seq = nt_block[target_raw_seq_idx]
                                nt_sequence_overall_weights[pfm_idx][nt_seq] = \
                                    nt_sequence_overall_weights[pfm_idx].get(nt_seq, 0.0) + norm_w
                # else: print("  Skipping NT PFM update (PFMs not init or nt_region_len <=0).")


                if aa_region_len_from_args > 0:
                    current_translated_aa_block_selections, det_aa_len_this_batch = translate_sequence_blocks(
                        current_raw_nt_blocks, selected_indices_0based_in_block,
                        aa_region_len_from_args, n_total_seqs_per_block
                    )
                    
                    if not aa_pfms_initialized_with_correct_length: # Try to initialize AA PFMs
                        if det_aa_len_this_batch is not None and det_aa_len_this_batch > 0:
                            actual_aa_len_for_pfm_established = det_aa_len_this_batch
                            print(f"  Established consistent AA length for PFM processing: {actual_aa_len_for_pfm_established}")
                            overall_aa_pfms = initialize_pfms(num_logos_to_process, amino_acids, actual_aa_len_for_pfm_established)
                            if overall_aa_pfms is None or (not overall_aa_pfms and num_logos_to_process > 0):
                                print("  Error: Failed to initialize AA PFMs. AA processing skipped.")
                                aa_pfms_initialized_with_correct_length = False; overall_aa_pfms = []
                            else:
                                aa_pfms_initialized_with_correct_length = True
                        # else: remain False, wait for a file that gives a valid det_aa_len_this_batch
                    
                    if aa_pfms_initialized_with_correct_length and overall_aa_pfms:
                        if det_aa_len_this_batch is not None and det_aa_len_this_batch != actual_aa_len_for_pfm_established:
                            print(f"  Warn: AA len in current batch ({det_aa_len_this_batch}) inconsistent with established ({actual_aa_len_for_pfm_established}) for '{fpath.name}'. Skipping AA PFM/freq for this file.")
                        elif not current_translated_aa_block_selections:
                             pass # print(f"  No valid AA blocks from '{fpath.name}' for PFM/freq.")
                        else:
                            update_pfms(overall_aa_pfms, current_translated_aa_block_selections, norm_w, amino_acids,
                                        actual_aa_len_for_pfm_established, is_aa_data=True)
                            for pfm_idx in range(num_logos_to_process):
                                for aa_selection_list_for_one_raw_block in current_translated_aa_block_selections:
                                    if pfm_idx < len(aa_selection_list_for_one_raw_block):
                                        aa_seq = aa_selection_list_for_one_raw_block[pfm_idx]
                                        if aa_seq is not None and len(aa_seq) == actual_aa_len_for_pfm_established:
                                            aa_sequence_overall_weights[pfm_idx][aa_seq] = \
                                                aa_sequence_overall_weights[pfm_idx].get(aa_seq, 0.0) + norm_w
                    elif not aa_pfms_initialized_with_correct_length and aa_region_len_from_args > 0 and file_idx == len(log_w_data) -1 :
                         print("  Warn: End of files, but could not establish consistent AA length. AA logos/freqs not generated.")
                elif file_idx == 0: print("  Skipping AA processing (aa_start/aa_end implies zero or negative length).")

        except FileNotFoundError: print(f"\nError: File disappeared during Pass 2: {fpath}")
        except Exception as e: print(f"\nError Pass 2 processing '{fpath.name}': {e}"); traceback.print_exc()

    print(f"\n--- Finished processing {files_proc_p2} files in Pass 2 ---")
    print(f"Sum of normalized weights processed: {total_norm_w_sum:.6f}")

    if nt_sequence_overall_weights and any(nt_sequence_overall_weights):
        # --- MAP (Highest Frequency) Report for NT ---
        find_and_print_highest_frequency_sequences(nt_sequence_overall_weights, "Nucleotide", plot_title_indices_1based)
        # --- NEW: MM (Marginal Mode) Report for NT ---
        find_and_print_marginal_mode_sequences(overall_nt_pfms, "Nucleotide", plot_title_indices_1based)
    else:
        print("\nNo Nucleotide sequence frequency data to report for MAP/MM.")

    if aa_pfms_initialized_with_correct_length and aa_sequence_overall_weights and any(aa_sequence_overall_weights):
        # --- MAP (Highest Frequency) Report for AA ---
        find_and_print_highest_frequency_sequences(aa_sequence_overall_weights, "Amino Acid", plot_title_indices_1based, start_aa_arg)
        # --- NEW: MM (Marginal Mode) Report for AA ---
        find_and_print_marginal_mode_sequences(overall_aa_pfms, "Amino Acid", plot_title_indices_1based, start_aa_arg)
    elif aa_region_len_from_args > 0:
        print("\nNo Amino Acid sequence frequency data (no consistent AA length or no valid sequences).")

    nt_plotted = False; aa_plotted = False
    base_nt_fname_str = ""; base_aa_fname_str = ""

    if overall_nt_pfms and any(not pfm.empty for pfm in overall_nt_pfms) and files_proc_p2 > 0 and nt_region_len > 0 :
        prob_mats_nt = generate_probability_matrices(overall_nt_pfms)
        if prob_mats_nt and any(not mat.empty for mat in prob_mats_nt):
            indices_str = "_".join(map(str,plot_title_indices_1based))
            base_nt_fname_str = f"{args.output_prefix}_nt_normW_AAslice{start_aa_arg}-{end_aa_arg}.pdf"
            plot_sequence_logos(prob_mats_nt, base_nt_fname_str,
                                plot_title_indices_1based=plot_title_indices_1based, is_aa=False)
            nt_plotted = True
        else: print("Failed to generate Nt probability matrices or all were empty.")
    else: print("No valid Nt PFM data. Nt logos not generated.")

    if aa_pfms_initialized_with_correct_length and overall_aa_pfms and any(not pfm.empty for pfm in overall_aa_pfms) and files_proc_p2 > 0 and actual_aa_len_for_pfm_established is not None and actual_aa_len_for_pfm_established > 0:
        prob_mats_aa = generate_probability_matrices(overall_aa_pfms)
        if prob_mats_aa and any(not mat.empty for mat in prob_mats_aa):
            indices_str = "_".join(map(str,plot_title_indices_1based))
            base_aa_fname_str = f"{args.output_prefix}_aa_normW_AAregion{start_aa_arg}-{end_aa_arg}.pdf"
            plot_sequence_logos(prob_mats_aa, base_aa_fname_str,
                                plot_title_indices_1based=plot_title_indices_1based,
                                is_aa=True, start_lbl_for_aa_axis=start_aa_arg)
            aa_plotted = True
        else: print("Failed to generate AA probability matrices or all were empty.")
    elif aa_region_len_from_args > 0:
        print("No valid AA PFM data. AA logos not generated.")

    print("\nProcessing finished.")
    if nt_plotted: print(f"Nt logo(s) pattern: '{Path(base_nt_fname_str).stem}_partX.pdf' or '{Path(base_nt_fname_str).name}'")
    if aa_plotted: print(f"AA logo(s) pattern: '{Path(base_aa_fname_str).stem}_partX.pdf' or '{Path(base_aa_fname_str).name}'")

if __name__ == "__main__":
    main()
