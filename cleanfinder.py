#!/usr/bin/env python3
"""
CleanFinder - CRISPR Editing Analysis Tool (CLI Version)
=========================================================
Analyzes CRISPR editing outcomes from FASTQ/FASTA files.

Developed by GEMD and EACR, IUF Leibniz Research Institute, Duesseldorf

Usage:
    python cleanfinder.py -r reference.fa -q reads.fastq -g GRNA_SEQUENCE

For help:
    python cleanfinder.py --help
"""

import argparse
import sys
import csv
import json
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path

# ============================================================================
# ASCII ART LOGO & WELCOME BANNER
# ============================================================================

VERSION = "1.0.0"

LOGO = r"""
   ____ _                  _____ _           _           
  / ___| | ___  __ _ _ __ |  ___(_)_ __   __| | ___ _ __ 
 | |   | |/ _ \/ _` | '_ \| |_  | | '_ \ / _` |/ _ \ '__|
 | |___| |  __/ (_| | | | |  _| | | | | | (_| |  __/ |   
  \____|_|\___|\__,_|_| |_|_|   |_|_| |_|\__,_|\___|_|   
                                                          
              CRISPR Editing Analysis Tool v{version}              
    
  MODES:    Cas9 | Cas12 | Prime Editing | Base Editing
  
  FEATURES:
    - Affine gap alignment (matches web app algorithm)
    - Fuzzy anchor matching with configurable mismatches
    - Phred quality filtering
    - In-frame / Out-of-frame detection
    - CSV export with aligned sequences
    - Configurable analysis window size
    
  Developed by GEMD and EACR, IUF Leibniz Research Institute, Duesseldorf
""".format(version=VERSION)

# ============================================================================
# SEQUENCE UTILITIES
# ============================================================================

def sanitize_seq(seq):
    """Remove non-DNA characters and convert to uppercase."""
    return ''.join(c for c in seq.upper() if c in 'ACGTN')


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(c, 'N') for c in reversed(seq))


# ============================================================================
# FILE PARSING
# ============================================================================

def parse_fasta(filepath):
    """Parse a FASTA file and return the first sequence."""
    seq_parts = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_parts:
                    break  # Only take first sequence
                continue
            seq_parts.append(line)
    return sanitize_seq(''.join(seq_parts))


def parse_fastq(filepath, min_phred=0):
    """Parse a FASTQ file and return sequences with quality filtering."""
    reads = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines) - 3:
        line = lines[i].strip()
        if line.startswith('@'):
            seq = sanitize_seq(lines[i + 1].strip())
            qual = lines[i + 3].strip()
            
            # Calculate average Phred score
            if qual and len(qual) > 0:
                avg_phred = sum(ord(c) - 33 for c in qual) / len(qual)
            else:
                avg_phred = 0
            
            # Apply quality filter
            if avg_phred >= min_phred and seq:
                reads.append({'seq': seq, 'phred': avg_phred})
            
            i += 4
        else:
            i += 1
    
    return reads


def parse_reference(filepath):
    """Parse reference from FASTA file or plain text."""
    path = Path(filepath)
    
    if not path.exists():
        # Maybe it's a direct sequence
        seq = sanitize_seq(filepath)
        if len(seq) > 10:
            return seq
        raise FileNotFoundError(f"Reference file not found: {filepath}")
    
    # Try FASTA format
    with open(filepath, 'r') as f:
        content = f.read()
    
    if content.startswith('>'):
        return parse_fasta(filepath)
    else:
        # Plain text
        return sanitize_seq(content)


# ============================================================================
# ENSEMBL REFERENCE FETCH
# ============================================================================

ENSEMBL_SPECIES = {
    'human':  'homo_sapiens',
    'mouse':  'mus_musculus',
    'rat':    'rattus_norvegicus',
    'zebrafish': 'danio_rerio',
    # Allow passing the full species name directly
}

def fetch_ensembl_reference(gene, grna=None, species='homo_sapiens', pad5=50, pad3=50, anchor_len=20):
    """
    Fetch a trimmed reference sequence from Ensembl REST API.
    Mirrors the web app Ensembl fetch logic exactly.

    Args:
        gene:       Gene symbol (e.g. 'SNAI1', 'BRCA1')
        grna:       gRNA sequence for smart trimming (optional but recommended)
        species:    Ensembl species string (default: homo_sapiens)
        pad5:       bp to include 5' of gRNA (default: 50)
        pad3:       bp to include 3' of gRNA (default: 50)
        anchor_len: extra bp on each side for anchors (default: 20)

    Returns:
        Trimmed reference sequence string (forward strand, uppercase)
    """
    # Resolve short species aliases
    species = ENSEMBL_SPECIES.get(species.lower(), species)

    BASE = 'https://rest.ensembl.org'
    HEADERS = {'Content-Type': 'application/json'}

    def api_get(url):
        req = urllib.request.Request(url, headers=HEADERS)
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            raise RuntimeError(f"Ensembl API error {e.code}: {url}")
        except urllib.error.URLError as e:
            raise RuntimeError(f"Network error: {e.reason}. Check your internet connection.")

    # Step 1 — Gene lookup
    print(f"  Querying Ensembl for '{gene}' ({species})...")
    lookup = api_get(f"{BASE}/lookup/symbol/{species}/{gene}?content-type=application/json")
    gene_id = lookup.get('id')
    if not gene_id:
        raise RuntimeError(f"Gene '{gene}' not found in Ensembl for species '{species}'.")
    print(f"  Found: {gene_id}  chr{lookup.get('seq_region_name')}:{lookup.get('start')}-{lookup.get('end')} ({lookup.get('strand',1)>0 and '+' or '-'})")

    # Step 2 — Sequence fetch with 1 kb flanks (same as web app)
    print(f"  Fetching sequence (±1kb flanks)...")
    seq_data = api_get(f"{BASE}/sequence/id/{gene_id}?content-type=application/json&expand_5prime=1000&expand_3prime=1000")
    full_seq = sanitize_seq(seq_data.get('seq', ''))
    if len(full_seq) < 10:
        raise RuntimeError(f"No sequence returned for {gene_id}.")
    print(f"  Raw sequence: {len(full_seq)} bp")

    # Step 3 — Smart trim around gRNA (mirrors web app logic)
    if grna:
        grna_clean = sanitize_seq(grna)
        match_idx = full_seq.find(grna_clean)
        match_rc  = False

        if match_idx == -1:
            rc = reverse_complement(grna_clean)
            match_idx = full_seq.find(rc)
            match_rc = True

        if match_idx == -1:
            print(f"  ⚠  gRNA not found in fetched sequence — returning full flanked gene.")
            return full_seq

        if match_rc:
            start = max(0, match_idx - pad3 - anchor_len)
            end   = min(len(full_seq), match_idx + len(grna_clean) + pad5 + anchor_len)
            trimmed = full_seq[start:end]
            trimmed = reverse_complement(trimmed)   # flip to forward
            strand_note = 'RC → flipped to forward'
        else:
            start = max(0, match_idx - pad5 - anchor_len)
            end   = min(len(full_seq), match_idx + len(grna_clean) + pad3 + anchor_len)
            trimmed = full_seq[start:end]
            strand_note = 'forward'

        print(f"  gRNA matched ({strand_note}). Trimmed to {len(trimmed)} bp.")
        return trimmed

    # No gRNA: return 2 kb window centred on the gene body
    centre  = len(full_seq) // 2
    half_kb = 1000
    trimmed = full_seq[max(0, centre - half_kb): centre + half_kb]
    print(f"  No gRNA provided — returning central 2 kb window ({len(trimmed)} bp).")
    return trimmed


# ============================================================================
# ANCHOR MATCHING
# ============================================================================

def fuzzy_index_of(haystack, needle, max_mismatches, start=0):
    """Find needle in haystack allowing up to max_mismatches mismatches."""
    if not needle:
        return -1
    
    h_len = len(haystack)
    n_len = len(needle)
    
    for i in range(start, h_len - n_len + 1):
        mismatches = 0
        for j in range(n_len):
            if haystack[i + j] != needle[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            return i
    
    return -1


def extract_middle(read, left_anchor, right_anchor, max_mismatches):
    """Extract the middle region between anchors."""
    # Try forward orientation
    l_pos = fuzzy_index_of(read, left_anchor, max_mismatches)
    if l_pos != -1:
        l_end = l_pos + len(left_anchor)
        r_pos = fuzzy_index_of(read, right_anchor, max_mismatches, l_end)
        if r_pos != -1:
            return read[l_end:r_pos], False  # middle, is_reverse
    
    # Try reverse complement
    left_rc = reverse_complement(right_anchor)
    right_rc = reverse_complement(left_anchor)
    
    l_pos = fuzzy_index_of(read, left_rc, max_mismatches)
    if l_pos != -1:
        l_end = l_pos + len(left_rc)
        r_pos = fuzzy_index_of(read, right_rc, max_mismatches, l_end)
        if r_pos != -1:
            middle_rc = read[l_end:r_pos]
            return reverse_complement(middle_rc), True  # Convert back to forward
    
    return None, None


# ============================================================================
# ALIGNMENT (Affine Gap Penalty - Matching Web App)
# ============================================================================

def fitting_align(ref, query, match_score=2, mismatch_score=-6, gap_open=-10, gap_ext=-2):
    """
    Perform fitting alignment using affine gap penalties.
    Uses 3-state DP (M/Ix/Iy) matching the web app's algorithm.
    Returns alignment strings and score.
    """
    m = len(ref)
    n = len(query)
    
    if m == 0 or n == 0:
        return {'alnRef': ref, 'alnRead': query, 'score': 0, 'ok': False}
    
    NEG_INF = -999999999
    
    # Initialize 3 DP matrices: M (match), Ix (gap in query), Iy (gap in ref)
    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Ix = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Iy = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    
    # Initialize (0,0)
    M[0][0] = 0
    
    # Initialize first row (free start in query for fitting alignment)
    for j in range(1, n + 1):
        M[0][j] = 0  # Free start in query
    
    # Initialize first column (gap penalties for ref)
    for i in range(1, m + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_ext
    
    # Fill DP tables
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Match/mismatch score
            score = match_score if ref[i-1] == query[j-1] else mismatch_score
            
            # Iy: gap in ref (insertion in read)
            Iy[i][j] = max(M[i][j-1] + gap_open, Iy[i][j-1] + gap_ext)
            
            # Ix: gap in query (deletion in read)
            Ix[i][j] = max(M[i-1][j] + gap_open, Ix[i-1][j] + gap_ext)
            
            # M: match/mismatch
            best_prev = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
            M[i][j] = best_prev + score
    
    # Find best ending position (fitting: ref fully consumed, query anywhere)
    best_score = NEG_INF
    best_j = n
    for j in range(n + 1):
        for mat in [M[m][j], Ix[m][j]]:
            if mat > best_score:
                best_score = mat
                best_j = j
    
    # Traceback
    aln_ref = []
    aln_read = []
    i, j = m, best_j
    
    # Add trailing unaligned query bases
    while j < n:
        aln_ref.append('-')
        aln_read.append(query[j])
        j += 1
    
    # Determine starting state
    if M[i][best_j] >= Ix[i][best_j]:
        state = 'M'
    else:
        state = 'Ix'
    
    j = best_j
    
    # Main traceback
    while i > 0 and j > 0:
        if state == 'M':
            score = match_score if ref[i-1] == query[j-1] else mismatch_score
            aln_ref.append(ref[i-1])
            aln_read.append(query[j-1])
            
            # Determine where we came from
            val = M[i][j] - score
            if M[i-1][j-1] == val:
                state = 'M'
            elif Ix[i-1][j-1] == val:
                state = 'Ix'
            else:
                state = 'Iy'
            i -= 1
            j -= 1
            
        elif state == 'Ix':
            aln_ref.append(ref[i-1])
            aln_read.append('-')
            
            # Determine where we came from
            if Ix[i][j] == Ix[i-1][j] + gap_ext:
                state = 'Ix'
            else:
                state = 'M'
            i -= 1
            
        else:  # state == 'Iy'
            aln_ref.append('-')
            aln_read.append(query[j-1])
            
            if Iy[i][j] == Iy[i][j-1] + gap_ext:
                state = 'Iy'
            else:
                state = 'M'
            j -= 1
    
    # Handle remaining characters
    while j > 0:
        aln_ref.append('-')
        aln_read.append(query[j-1])
        j -= 1
    
    while i > 0:
        aln_ref.append(ref[i-1])
        aln_read.append('-')
        i -= 1
    
    # Reverse (we built backwards)
    aln_ref_str = ''.join(reversed(aln_ref))
    aln_read_str = ''.join(reversed(aln_read))
    
    return {
        'alnRef': aln_ref_str,
        'alnRead': aln_read_str,
        'score': best_score,
        'ok': True
    }


def detect_microhomology(aln_ref, aln_read):
    """
    Detect microhomology at deletion junctions.
    Returns: {'length': int, 'sequence': str, 'pathway': 'MMEJ'|'NHEJ'}
    """
    # Find the deletion region (consecutive gaps in read)
    gap_start = -1
    gap_end = -1
    
    for i, (r, q) in enumerate(zip(aln_ref, aln_read)):
        if q == '-':
            if gap_start == -1:
                gap_start = i
            gap_end = i
    
    if gap_start == -1:
        return None  # No deletion
    
    # Get the deleted sequence
    deleted_seq = aln_ref[gap_start:gap_end + 1].replace('-', '')
    
    if len(deleted_seq) < 2:
        return {'length': 0, 'sequence': '', 'pathway': 'NHEJ'}
    
    # Check for microhomology at junction
    # Look for sequence before deletion that matches sequence after deletion
    pre_seq = aln_ref[max(0, gap_start - 10):gap_start].replace('-', '')
    post_seq = aln_ref[gap_end + 1:gap_end + 11].replace('-', '')
    
    mh_length = 0
    mh_seq = ''
    
    # Check if bases at deletion end match bases after deletion (forward MH)
    # This is the classic microhomology pattern
    for length in range(1, min(len(deleted_seq), 10) + 1):
        # Check if end of deleted sequence matches start of post-sequence
        if length <= len(post_seq) and deleted_seq[-length:] == post_seq[:length]:
            if length > mh_length:
                mh_length = length
                mh_seq = deleted_seq[-length:]
        # Check if start of deleted sequence matches end of pre-sequence
        if length <= len(pre_seq) and deleted_seq[:length] == pre_seq[-length:]:
            if length > mh_length:
                mh_length = length
                mh_seq = deleted_seq[:length]
    
    pathway = 'MMEJ' if mh_length >= 2 else 'NHEJ'
    
    return {
        'length': mh_length,
        'sequence': mh_seq,
        'pathway': pathway,
        'del_start': gap_start,
        'del_end': gap_end
    }


# ============================================================================
# CLASSIFICATION
# ============================================================================

def classify_genotype(middle, ref_middle, mode='cas9', hdr_template=None, 
                      pe_edit=None, be_conversion=None):
    """Classify a genotype based on the extracted middle sequence."""
    delta = len(middle) - len(ref_middle)
    
    # Check for Wildtype
    if middle == ref_middle:
        return {'type': 'WT', 'delta': 0, 'frame': '-'}
    
    # Check for Prime Edit (if PE mode and pe_edit provided)
    if mode == 'prime' and pe_edit:
        if pe_edit in middle or reverse_complement(pe_edit) in middle:
            return {
                'type': 'PE',
                'delta': delta,
                'frame': 'In-Frame' if delta % 3 == 0 else 'Out-Frame'
            }
    
    # Check for HDR/Knock-in
    if hdr_template and hdr_template in middle:
        return {
            'type': 'KI',
            'delta': delta,
            'frame': 'In-Frame' if delta % 3 == 0 else 'Out-Frame'
        }
    
    # Check for Base Edit (if BE mode)
    if mode == 'base' and be_conversion and delta == 0:
        # Count expected conversions
        from_base, to_base = be_conversion.upper().split('>') if '>' in be_conversion else ('C', 'T')
        conversions = 0
        for r, m in zip(ref_middle, middle):
            if r == from_base and m == to_base:
                conversions += 1
        if conversions > 0:
            return {'type': 'BE', 'delta': 0, 'frame': '-', 'conversions': conversions}
    
    # Standard classification
    if delta == 0:
        return {'type': 'SNP', 'delta': 0, 'frame': '-'}
    elif delta < 0:
        frame = 'In-Frame' if delta % 3 == 0 else 'Out-Frame'
        return {'type': 'DEL', 'delta': delta, 'frame': frame}
    else:
        frame = 'In-Frame' if delta % 3 == 0 else 'Out-Frame'
        return {'type': 'INS', 'delta': delta, 'frame': frame}


# ============================================================================
# ANALYSIS
# ============================================================================

def get_alignment_params(mode):
    """Get alignment parameters based on analysis mode."""
    params = {
        'cas9': {'match': 2, 'mismatch': -6, 'gap_open': -10, 'gap_ext': -2},
        'cas12': {'match': 2, 'mismatch': -6, 'gap_open': -10, 'gap_ext': -2},
        'prime': {'match': 2, 'mismatch': -3, 'gap_open': -1, 'gap_ext': -0.5},
        'base': {'match': 2, 'mismatch': -6, 'gap_open': -20, 'gap_ext': -20}
    }
    return params.get(mode, params['cas9'])


def analyze_reads(reads, ref_middle, left_anchor, right_anchor, 
                  max_mismatches=2, mode='cas9', hdr_template=None,
                  pe_edit=None, be_conversion=None):
    """Analyze all reads and return genotype counts."""
    genotypes = defaultdict(lambda: {
        'count': 0,
        'middle': '',
        'classification': None,
        'alignment': None
    })
    
    aligned = 0
    unaligned = 0
    
    # Get mode-specific alignment parameters
    align_params = get_alignment_params(mode)
    
    for read_data in reads:
        seq = read_data['seq']
        
        # Extract middle region
        middle, is_reverse = extract_middle(seq, left_anchor, right_anchor, max_mismatches)
        
        if middle is None:
            unaligned += 1
            continue
        
        aligned += 1
        
        # Use middle sequence as signature
        sig = middle
        
        if genotypes[sig]['count'] == 0:
            # First time seeing this genotype
            genotypes[sig]['middle'] = middle
            genotypes[sig]['classification'] = classify_genotype(
                middle, ref_middle, mode, hdr_template, pe_edit, be_conversion
            )
            genotypes[sig]['alignment'] = fitting_align(
                ref_middle, middle,
                align_params['match'],
                align_params['mismatch'],
                int(align_params['gap_open']),
                int(align_params['gap_ext'])
            )
            
            # Detect microhomology for deletions
            if genotypes[sig]['classification']['type'] == 'DEL' and genotypes[sig]['alignment']['ok']:
                mh_info = detect_microhomology(
                    genotypes[sig]['alignment']['alnRef'],
                    genotypes[sig]['alignment']['alnRead']
                )
                genotypes[sig]['microhomology'] = mh_info
            else:
                genotypes[sig]['microhomology'] = None
        
        genotypes[sig]['count'] += 1
    
    return {
        'genotypes': dict(genotypes),
        'aligned': aligned,
        'unaligned': unaligned,
        'total': len(reads)
    }


# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

def format_alignment_display(ref_aln, read_aln, grna=None, width=120):
    """Format alignment for text display with markers."""
    lines = []
    
    # Build marker line
    markers = []
    for r, q in zip(ref_aln, read_aln):
        if r == q:
            markers.append(' ')
        elif r == '-':
            markers.append('+')  # Insertion
        elif q == '-':
            markers.append('^')  # Deletion
        else:
            markers.append('*')  # Mismatch
    
    marker_str = ''.join(markers)
    
    # Split into chunks for display
    for i in range(0, len(ref_aln), width):
        chunk_ref = ref_aln[i:i+width]
        chunk_read = read_aln[i:i+width]
        chunk_mark = marker_str[i:i+width]
        
        lines.append(f"Ref:  {chunk_ref}")
        lines.append(f"Read: {chunk_read}")
        if '+' in chunk_mark or '^' in chunk_mark or '*' in chunk_mark:
            lines.append(f"      {chunk_mark}")
        lines.append("")
    
    return '\n'.join(lines)


def print_results(results, ref_middle, mode, show_alignment=True, top_n=10):
    """Print analysis results to console."""
    genotypes = results['genotypes']
    
    # Sort by count
    sorted_genos = sorted(genotypes.items(), key=lambda x: x[1]['count'], reverse=True)
    
    # Calculate stats
    total_aligned = results['aligned']
    wt_count = sum(g['count'] for _, g in sorted_genos if g['classification']['type'] == 'WT')
    edited_count = total_aligned - wt_count
    
    # Print header
    print(LOGO)
    print("=" * 70)
    print(f"Mode: {mode.upper()}")
    print("=" * 70)
    print(f"\nReference: {len(ref_middle)} bp (middle region)")
    print(f"Total reads: {results['total']:,}")
    print(f"Aligned: {results['aligned']:,} ({100*results['aligned']/max(results['total'],1):.1f}%)")
    print(f"Unaligned: {results['unaligned']:,}")
    
    print(f"\n--- Summary ---")
    print(f"  Wildtype:  {wt_count:,} ({100*wt_count/max(total_aligned,1):.1f}%)")
    print(f"  Edited:    {edited_count:,} ({100*edited_count/max(total_aligned,1):.1f}%)")
    
    # Count by type
    type_counts = defaultdict(int)
    frame_counts = {'In-Frame': 0, 'Out-Frame': 0}
    for _, g in sorted_genos:
        c = g['classification']
        type_counts[c['type']] += g['count']
        if c['frame'] in frame_counts:
            frame_counts[c['frame']] += g['count']
    
    # Helper function for bar charts
    def make_bar(value, max_value, width=30):
        if max_value == 0:
            return '░' * width
        filled = int((value / max_value) * width)
        return '█' * filled + '░' * (width - filled)
    
    if edited_count > 0:
        print(f"\n{'═' * 50}")
        print("GENOTYPE DISTRIBUTION")
        print('═' * 50)
        
        # Show WT first, then others
        all_types = [('WT', type_counts['WT'])] + [(t, type_counts[t]) for t in ['DEL', 'INS', 'SNP', 'KI', 'PE', 'BE'] if type_counts[t] > 0]
        max_count = max(c for _, c in all_types) if all_types else 1
        
        for typ, count in all_types:
            pct = 100 * count / max(total_aligned, 1)
            bar = make_bar(count, max_count)
            print(f"  {typ:6} {bar} {pct:5.1f}%  ({count:,})")
        
        print(f"\n{'═' * 50}")
        print("FRAME STATUS (edits only)")
        print('═' * 50)
        
        in_pct = 100 * frame_counts['In-Frame'] / max(edited_count, 1)
        out_pct = 100 * frame_counts['Out-Frame'] / max(edited_count, 1)
        
        in_bar = make_bar(frame_counts['In-Frame'], edited_count)
        out_bar = make_bar(frame_counts['Out-Frame'], edited_count)
        
        print(f"  In-Frame  {in_bar} {in_pct:5.1f}%")
        print(f"  Out-Frame {out_bar} {out_pct:5.1f}%")
        
        # Count repair pathways for deletions
        mmej_count = 0
        nhej_count = 0
        for _, g in sorted_genos:
            if g['classification']['type'] == 'DEL' and g.get('microhomology'):
                if g['microhomology']['pathway'] == 'MMEJ':
                    mmej_count += g['count']
                else:
                    nhej_count += g['count']
        
        if mmej_count > 0 or nhej_count > 0:
            print(f"\n{'═' * 50}")
            print("REPAIR PATHWAYS (deletions only)")
            print('═' * 50)
            
            del_total = mmej_count + nhej_count
            mmej_pct = 100 * mmej_count / max(del_total, 1)
            nhej_pct = 100 * nhej_count / max(del_total, 1)
            
            mmej_bar = make_bar(mmej_count, del_total)
            nhej_bar = make_bar(nhej_count, del_total)
            
            print(f"  MMEJ  {mmej_bar} {mmej_pct:5.1f}%  (MH ≥2bp)")
            print(f"  NHEJ  {nhej_bar} {nhej_pct:5.1f}%  (MH <2bp)")
    
    # Print top genotypes with alignment
    print(f"\n" + "=" * 70)
    print(f"Top {min(top_n, len(sorted_genos))} Genotypes")
    print("=" * 70)
    
    for i, (sig, geno) in enumerate(sorted_genos[:top_n]):
        c = geno['classification']
        count = geno['count']
        pct = 100 * count / max(total_aligned, 1)
        
        delta_str = f"{c['delta']:+d}" if c['delta'] != 0 else "0"
        
        # Add repair pathway info for deletions
        pathway_str = ""
        if c['type'] == 'DEL' and geno.get('microhomology'):
            mh = geno['microhomology']
            pathway_str = f"  [{mh['pathway']}"
            if mh['length'] > 0:
                pathway_str += f" MH:{mh['length']}bp '{mh['sequence']}'"
            pathway_str += "]"
        
        print(f"\n#{i+1}  {c['type']} {delta_str} ({c['frame']}){pathway_str}  —  {count:,} reads ({pct:.1f}%)")
        print("-" * 70)
        
        if show_alignment and geno['alignment'] and geno['alignment']['ok']:
            aln = geno['alignment']
            print(format_alignment_display(aln['alnRef'], aln['alnRead']))
    
    print("\n" + "=" * 70)


def export_csv(results, output_path, ref_middle):
    """Export results to CSV file."""
    genotypes = results['genotypes']
    sorted_genos = sorted(genotypes.items(), key=lambda x: x[1]['count'], reverse=True)
    
    total_aligned = results['aligned']
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Type', 'Delta', 'Frame', 'Sequence', 'Aligned_Ref', 'Aligned_Read', 'Reads', 'Percent'])
        
        for sig, geno in sorted_genos:
            c = geno['classification']
            count = geno['count']
            pct = f"{100*count/max(total_aligned,1):.2f}%"
            
            aln = geno['alignment']
            aln_ref = aln['alnRef'] if aln and aln['ok'] else ''
            aln_read = aln['alnRead'] if aln and aln['ok'] else ''
            
            writer.writerow([
                c['type'],
                c['delta'],
                c['frame'],
                geno['middle'],
                aln_ref,
                aln_read,
                count,
                pct
            ])
    
    print(f"\nResults exported to: {output_path}")


def print_multiwell_summary(all_results):
    """Print a multiwell-style summary table for multiple samples."""
    print(f"\n{'═' * 80}")
    print("MULTIWELL SUMMARY")
    print('═' * 80)
    
    # Header
    print(f"{'Sample':<30} {'Reads':>10} {'Aligned':>10} {'Edited %':>10} {'In-Frame %':>12}")
    print('─' * 80)
    
    # Data rows
    for sample_name, result in all_results:
        total = result['total']
        aligned = result['aligned']
        
        # Calculate editing efficiency
        genotypes = result['genotypes']
        wt_count = sum(g['count'] for _, g in genotypes.items() if g['classification']['type'] == 'WT')
        edited = aligned - wt_count
        edited_pct = 100 * edited / max(aligned, 1)
        
        # Calculate in-frame
        in_frame = sum(g['count'] for _, g in genotypes.items() 
                      if g['classification']['frame'] == 'In-Frame')
        in_frame_pct = 100 * in_frame / max(edited, 1)
        
        # Truncate long sample names
        display_name = sample_name[:28] + '..' if len(sample_name) > 30 else sample_name
        
        print(f"{display_name:<30} {total:>10,} {aligned:>10,} {edited_pct:>9.1f}% {in_frame_pct:>11.1f}%")
    
    print('═' * 80)


def export_multiwell_csv(all_results, output_path):
    """Export multiwell summary to CSV."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Sample', 'Total_Reads', 'Aligned', 'Aligned_Pct', 'WT', 'Edited', 
                        'Edited_Pct', 'DEL', 'INS', 'SNP', 'In_Frame', 'Out_Frame', 
                        'MMEJ', 'NHEJ'])
        
        for sample_name, result in all_results:
            total = result['total']
            aligned = result['aligned']
            genotypes = result['genotypes']
            
            # Count by type
            wt = sum(g['count'] for _, g in genotypes.items() if g['classification']['type'] == 'WT')
            dels = sum(g['count'] for _, g in genotypes.items() if g['classification']['type'] == 'DEL')
            ins = sum(g['count'] for _, g in genotypes.items() if g['classification']['type'] == 'INS')
            snps = sum(g['count'] for _, g in genotypes.items() if g['classification']['type'] == 'SNP')
            
            edited = aligned - wt
            in_frame = sum(g['count'] for _, g in genotypes.items() 
                          if g['classification']['frame'] == 'In-Frame')
            out_frame = sum(g['count'] for _, g in genotypes.items() 
                           if g['classification']['frame'] == 'Out-Frame')
            
            # MMEJ/NHEJ counts
            mmej = sum(g['count'] for _, g in genotypes.items() 
                      if g.get('microhomology') and g['microhomology']['pathway'] == 'MMEJ')
            nhej = sum(g['count'] for _, g in genotypes.items() 
                      if g.get('microhomology') and g['microhomology']['pathway'] == 'NHEJ')
            
            aligned_pct = 100 * aligned / max(total, 1)
            edited_pct = 100 * edited / max(aligned, 1)
            
            writer.writerow([sample_name, total, aligned, f"{aligned_pct:.1f}", wt, edited,
                           f"{edited_pct:.1f}", dels, ins, snps, in_frame, out_frame,
                           mmej, nhej])
    
    print(f"\nMultiwell summary exported to: {output_path}")


def parse_batch_config(config_path):
    """Parse batch config CSV file."""
    samples = []
    with open(config_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append({
                'sample_name': row.get('sample_name', row.get('name', '')),
                'fastq': row.get('fastq', row.get('reads', '')),
                'reference': row.get('reference', row.get('ref', '')),
                'grna': row.get('grna', row.get('gRNA', '')),
                'mode': row.get('mode', 'cas9').lower()
            })
    return samples


def analyze_single_file(fastq_path, ref_seq, grna, args, sample_name=None):
    """Analyze a single FASTQ file and return results."""
    import os
    if sample_name is None:
        sample_name = os.path.basename(fastq_path)
    
    anchor_len = args.anchor_len
    
    # Define anchors based on gRNA
    if grna:
        grna = sanitize_seq(grna)
        grna_pos = ref_seq.find(grna)
        if grna_pos == -1:
            grna_pos = ref_seq.find(reverse_complement(grna))
        
        if grna_pos != -1:
            center = grna_pos + len(grna) // 2
            half_window = args.window_size
            
            max_half_window = min(center, len(ref_seq) - center) - anchor_len
            if half_window > max_half_window:
                half_window = max(20, max_half_window)
            
            window_start = max(0, center - half_window)
            window_end = min(len(ref_seq), center + half_window)
            
            left_anchor = ref_seq[max(0, window_start - anchor_len):window_start]
            right_anchor = ref_seq[window_end:min(len(ref_seq), window_end + anchor_len)]
            ref_middle = ref_seq[window_start:window_end]
        else:
            left_anchor = ref_seq[:anchor_len]
            right_anchor = ref_seq[-anchor_len:]
            ref_middle = ref_seq[anchor_len:-anchor_len]
    else:
        left_anchor = ref_seq[:anchor_len]
        right_anchor = ref_seq[-anchor_len:]
        ref_middle = ref_seq[anchor_len:-anchor_len]
    
    # Parse reads
    reads = parse_fastq(fastq_path, args.min_phred)
    
    # Analyze
    hdr_template = sanitize_seq(args.hdr) if args.hdr else None
    pe_edit = sanitize_seq(args.pe_edit) if args.pe_edit else None
    
    results = analyze_reads(
        reads, ref_middle, left_anchor, right_anchor,
        max_mismatches=args.max_mismatches,
        mode=args.mode,
        hdr_template=hdr_template,
        pe_edit=pe_edit,
        be_conversion=args.be_conversion
    )
    
    return results, ref_middle


# ============================================================================
# MAIN
# ============================================================================

def main():
    # Show welcome banner when running without arguments
    if len(sys.argv) == 1:
        print(LOGO)
        print("Usage: python cleanfinder.py -r REFERENCE -q READS.fastq [options]")
        print("\nRun 'python cleanfinder.py --help' for full options.\n")
        sys.exit(0)
    
    parser = argparse.ArgumentParser(
        description='CleanFinder - CRISPR Editing Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard Cas9 analysis
  python cleanfinder.py -r reference.fa -q reads.fastq -g ATCGATCGATCG

  # Prime Editing mode
  python cleanfinder.py -r ref.fa -q reads.fq --mode prime --pe-edit EXPECTED_EDIT

  # Base Editing mode
  python cleanfinder.py -r ref.fa -q reads.fq --mode base --be-conversion C>T

  # With quality filtering
  python cleanfinder.py -r ref.fa -q sample.fq -o results.csv --min-phred 20

  # Multi-file analysis (same reference, multiple samples)
  python cleanfinder.py -r ref.fa -q sample1.fq sample2.fq sample3.fq -g GRNA

  # Batch analysis from config file
  python cleanfinder.py --batch config.csv -o results/
  
  Config CSV format: sample_name,fastq,reference,grna,mode

Developed by GEMD and EACR, IUF Leibniz Research Institute, Duesseldorf
        """
    )
    
    # Main arguments
    parser.add_argument('-r', '--ref', default=None,
                        help='Reference sequence (FASTA file or plain text). Overridden by --gene.')
    parser.add_argument('-q', '--reads', nargs='+', default=None,
                        help='FASTQ file(s) with sequencing reads (supports multiple files)')
    parser.add_argument('-g', '--grna', default=None,
                        help='gRNA sequence for window anchoring (optional)')
    parser.add_argument('-o', '--output', default='results.csv',
                        help='Output CSV file or directory for batch mode (default: results.csv)')

    # Ensembl fetch
    parser.add_argument('--gene', default=None,
                        help='Fetch reference from Ensembl by gene symbol (e.g. SNAI1). Requires internet. Overrides -r.')
    parser.add_argument('--species', default='homo_sapiens',
                        help='Ensembl species for --gene (default: homo_sapiens). Aliases: human, mouse, rat, zebrafish.')
    
    # Batch mode
    parser.add_argument('--batch', default=None,
                        help='CSV config file for batch analysis (columns: sample_name,fastq,reference,grna,mode)')
    
    # Mode settings
    parser.add_argument('--mode', choices=['cas9', 'cas12', 'prime', 'base'],
                        default='cas9', help='Analysis mode (default: cas9)')
    parser.add_argument('--hdr', default=None,
                        help='HDR/Knock-in template sequence (optional)')
    parser.add_argument('--pe-edit', default=None,
                        help='Expected Prime Edit sequence (for PE mode)')
    parser.add_argument('--be-conversion', default=None,
                        help='Base conversion pattern e.g. C>T (for BE mode)')
    
    # Alignment parameters
    parser.add_argument('--anchor-len', type=int, default=20,
                        help='Anchor length in bp (default: 20)')
    parser.add_argument('--max-mismatches', type=int, default=2,
                        help='Max mismatches for anchor matching (default: 2)')
    parser.add_argument('--min-phred', type=int, default=0,
                        help='Minimum Phred score filter (default: 0)')
    parser.add_argument('--window-size', type=int, default=40,
                        help='Window size on each side of gRNA in bp (default: 40, max ~200)')
    
    # Output options
    parser.add_argument('--top', type=int, default=15,
                        help='Number of top genotypes to display (default: 15)')
    parser.add_argument('--show-all', action='store_true',
                        help='Show all genotypes in terminal and PDF (overrides --top)')
    parser.add_argument('--no-alignment', action='store_true',
                        help='Skip alignment display in output')
    parser.add_argument('--quiet', action='store_true',
                        help='Minimal output (useful for batch mode)')
    parser.add_argument('--pdf', action='store_true',
                        help='Generate a PDF report with charts and alignment figures (requires matplotlib)')
    
    args = parser.parse_args()
    
    # =========================================================================
    # BATCH MODE: Process multiple samples from config CSV
    # =========================================================================
    if args.batch:
        import os
        print(LOGO)
        print(f"═══════════════════════════════════════════════════════════════════════")
        print(f"BATCH MODE - Processing config: {args.batch}")
        print(f"═══════════════════════════════════════════════════════════════════════\n")
        
        try:
            samples = parse_batch_config(args.batch)
        except Exception as e:
            print(f"Error reading batch config: {e}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(samples)} samples to process\n")
        
        # Create output directory if needed
        output_dir = args.output if args.output.endswith('/') else args.output + '/'
        os.makedirs(output_dir, exist_ok=True)
        
        all_results = []
        
        for i, sample in enumerate(samples):
            sample_name = sample['sample_name'] or os.path.basename(sample['fastq'])
            print(f"[{i+1}/{len(samples)}] Processing {sample_name}...")
            
            try:
                # Parse reference for this sample
                ref_seq = parse_reference(sample['reference'])
                
                # Override mode if specified in config
                args.mode = sample['mode'] or 'cas9'
                
                results, ref_middle = analyze_single_file(
                    sample['fastq'], ref_seq, sample['grna'], args, sample_name
                )
                
                all_results.append((sample_name, results))
                
                # Export individual CSV
                csv_path = f"{output_dir}{sample_name.replace('.fastq', '').replace('.fq', '')}.csv"
                export_csv(results, csv_path, ref_middle)
                
                # Print brief status
                wt = sum(g['count'] for _, g in results['genotypes'].items() 
                        if g['classification']['type'] == 'WT')
                edited_pct = 100 * (results['aligned'] - wt) / max(results['aligned'], 1)
                print(f"    ✓ {results['aligned']:,} aligned, {edited_pct:.1f}% edited")
                
            except Exception as e:
                print(f"    ✗ Error: {e}")
                continue
        
        # Print and export multiwell summary
        if all_results:
            print_multiwell_summary(all_results)
            summary_path = f"{output_dir}batch_summary.csv"
            export_multiwell_csv(all_results, summary_path)

            # Export per-sample PDF reports (optional)
            if args.pdf:
                for s_name, s_results in all_results:
                    top_sig = max(s_results['genotypes'], key=lambda k: s_results['genotypes'][k]['count'], default=None)
                    if top_sig:
                        aln = s_results['genotypes'][top_sig].get('alignment', {})
                        ref_m = aln.get('alnRef', '').replace('-', '') if aln and aln.get('ok') else ''
                    else:
                        ref_m = ''
                    pdf_path = f"{output_dir}{s_name.replace('.fastq','').replace('.fq','')}.pdf"
                    generate_pdf_report(s_results, ref_m, args.mode, s_name, pdf_path, args.show_all)

        print(f"\nBatch processing complete! Results saved to: {output_dir}")
        sys.exit(0)
    
    # =========================================================================
    # MULTI-FILE MODE: Same reference, multiple FASTQ files
    # =========================================================================
    if args.reads and len(args.reads) > 1:
        import os
        
        if not args.ref:
            print("Error: Reference (-r) is required for multi-file mode", file=sys.stderr)
            sys.exit(1)
        
        print(LOGO)
        print(f"═══════════════════════════════════════════════════════════════════════")
        print(f"MULTI-FILE MODE - Processing {len(args.reads)} samples")
        print(f"═══════════════════════════════════════════════════════════════════════\n")
        
        # Parse reference once
        try:
            ref_seq = parse_reference(args.ref)
        except Exception as e:
            print(f"Error loading reference: {e}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Reference: {len(ref_seq)} bp")
        if args.grna:
            print(f"gRNA: {args.grna}")
        print()
        
        all_results = []
        
        for i, fastq_path in enumerate(args.reads):
            sample_name = os.path.basename(fastq_path)
            print(f"[{i+1}/{len(args.reads)}] Processing {sample_name}...")
            
            try:
                results, ref_middle = analyze_single_file(
                    fastq_path, ref_seq, args.grna, args, sample_name
                )
                
                all_results.append((sample_name, results))
                
                # Brief status
                wt = sum(g['count'] for _, g in results['genotypes'].items() 
                        if g['classification']['type'] == 'WT')
                edited_pct = 100 * (results['aligned'] - wt) / max(results['aligned'], 1)
                print(f"    ✓ {results['aligned']:,} aligned, {edited_pct:.1f}% edited")
                
            except Exception as e:
                print(f"    ✗ Error: {e}")
                continue
        
        # Print and export multiwell summary
        if all_results:
            print_multiwell_summary(all_results)

            # Export multiwell CSV
            summary_path = args.output.replace('.csv', '_multiwell.csv') if args.output.endswith('.csv') else 'multiwell_summary.csv'
            export_multiwell_csv(all_results, summary_path)

            # Export per-sample PDF reports (optional)
            if args.pdf:
                import os
                for s_name, s_results in all_results:
                    # We need ref_middle per sample — re-derive from the last results
                    # (stored indirectly; use the aligned ref from top genotype as proxy)
                    top_sig = max(s_results['genotypes'], key=lambda k: s_results['genotypes'][k]['count'], default=None)
                    if top_sig:
                        aln = s_results['genotypes'][top_sig].get('alignment', {})
                        ref_m = aln.get('alnRef', '').replace('-', '') if aln and aln.get('ok') else ''
                    else:
                        ref_m = ''
                    pdf_path = os.path.join(
                        args.output if not args.output.endswith('.csv') else os.path.dirname(args.output) or '.',
                        f"{s_name.replace('.fastq','').replace('.fq','')}.pdf"
                    )
                    generate_pdf_report(s_results, ref_m, args.mode, s_name, pdf_path, args.show_all)

        print("\nDone!")
        sys.exit(0)
    
    # =========================================================================
    # SINGLE FILE MODE: Standard analysis
    # =========================================================================
    if not args.ref or not args.reads:
        print("Error: Both reference (-r) and reads (-q) are required", file=sys.stderr)
        print("Run 'python cleanfinder.py --help' for usage.")
        sys.exit(1)
    
    # Parse reference (Ensembl fetch takes priority over -r)
    print("Loading reference...")
    try:
        if args.gene:
            ref_seq = fetch_ensembl_reference(
                args.gene, grna=args.grna,
                species=args.species,
                pad5=args.window_size, pad3=args.window_size,
                anchor_len=args.anchor_len
            )
        elif args.ref:
            ref_seq = parse_reference(args.ref)
        else:
            print("Error: provide either -r/--ref or --gene", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"Error loading reference: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Reference: {len(ref_seq)} bp")
    
    # Define anchors and middle region
    anchor_len = args.anchor_len
    
    if args.grna:
        # Find gRNA position and center window around it
        grna = sanitize_seq(args.grna)
        grna_pos = ref_seq.find(grna)
        if grna_pos == -1:
            grna_pos = ref_seq.find(reverse_complement(grna))
        
        if grna_pos != -1:
            # Center on gRNA with user-defined window size
            center = grna_pos + len(grna) // 2
            half_window = args.window_size
            
            # AUTO-SCALE: Check if window size is too large for reference
            max_half_window = min(center, len(ref_seq) - center) - args.anchor_len
            if half_window > max_half_window:
                old_window = half_window
                half_window = max(20, max_half_window)  # Minimum 20bp
                print(f"  ⚠️  Window size auto-scaled from {old_window} to {half_window} (reference too short)")
            
            window_start = max(0, center - half_window)
            window_end = min(len(ref_seq), center + half_window)
            
            left_anchor = ref_seq[max(0, window_start - anchor_len):window_start]
            right_anchor = ref_seq[window_end:min(len(ref_seq), window_end + anchor_len)]
            ref_middle = ref_seq[window_start:window_end]
            
            print(f"  gRNA found at position {grna_pos}")
        else:
            print(f"  Warning: gRNA not found in reference, using ends as anchors")
            left_anchor = ref_seq[:anchor_len]
            right_anchor = ref_seq[-anchor_len:]
            ref_middle = ref_seq[anchor_len:-anchor_len]
    else:
        # Use ends as anchors
        left_anchor = ref_seq[:anchor_len]
        right_anchor = ref_seq[-anchor_len:]
        ref_middle = ref_seq[anchor_len:-anchor_len]
    
    print(f"  Left anchor: {left_anchor[:20]}{'...' if len(left_anchor) > 20 else ''}")
    print(f"  Right anchor: {right_anchor[:20]}{'...' if len(right_anchor) > 20 else ''}")
    print(f"  Middle region: {len(ref_middle)} bp")
    
    # Parse reads (single file)
    fastq_path = args.reads[0] if isinstance(args.reads, list) else args.reads
    print(f"\nLoading reads from {fastq_path}...")
    try:
        reads = parse_fastq(fastq_path, args.min_phred)
    except Exception as e:
        print(f"Error loading reads: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Loaded {len(reads):,} reads (after Phred filter >= {args.min_phred})")
    
    # Analyze
    print(f"\nAnalyzing reads... (Mode: {args.mode.upper()})")
    hdr_template = sanitize_seq(args.hdr) if args.hdr else None
    pe_edit = sanitize_seq(args.pe_edit) if args.pe_edit else None
    
    results = analyze_reads(
        reads, 
        ref_middle, 
        left_anchor, 
        right_anchor,
        max_mismatches=args.max_mismatches,
        mode=args.mode,
        hdr_template=hdr_template,
        pe_edit=pe_edit,
        be_conversion=args.be_conversion
    )
    
    # Print results
    print_results(results, ref_middle, args.mode,
                  show_alignment=not args.no_alignment,
                  top_n=args.top)
    
    # Export CSV
    export_csv(results, args.output, ref_middle)

    # Export PDF report (optional)
    if args.pdf:
        sample_name = Path(fastq_path).stem
        pdf_path = args.output.replace('.csv', '.pdf') if args.output.endswith('.csv') else args.output + '.pdf'
        generate_pdf_report(results, ref_middle, args.mode, sample_name, pdf_path, args.show_all)

    print("\nDone!")


# ============================================================================
# PDF REPORT GENERATION (requires matplotlib - optional)
# ============================================================================

def _check_matplotlib():
    """Check if matplotlib is available."""
    try:
        import matplotlib
        return True
    except ImportError:
        return False


def generate_pdf_report(results, ref_middle, mode, sample_name, output_path, show_all=False):
    """
    Generate a publication-quality PDF report with charts and alignment figures.
    Requires matplotlib. Gracefully skipped if not installed.

    Parameters:
        show_all: If True, show all genotypes; if False, show top 15
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.gridspec as gridspec
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.patches import FancyBboxPatch
        import datetime
    except ImportError:
        print("\n⚠  PDF report skipped: matplotlib not installed.")
        print("   Install it with:  pip install matplotlib")
        return

    genotypes  = results['genotypes']
    sorted_genos = sorted(genotypes.items(), key=lambda x: x[1]['count'], reverse=True)
    total      = results['total']
    aligned    = results['aligned']
    unaligned  = results['unaligned']

    # ── Aggregate counts ────────────────────────────────────────────────────
    from collections import defaultdict
    type_counts  = defaultdict(int)
    frame_counts = {'In-Frame': 0, 'Out-Frame': 0}
    mmej_count   = 0
    nhej_count   = 0

    for _, g in sorted_genos:
        c = g['classification']
        type_counts[c['type']] += g['count']
        if c['frame'] in frame_counts:
            frame_counts[c['frame']] += g['count']
        if g.get('microhomology'):
            if g['microhomology']['pathway'] == 'MMEJ':
                mmej_count += g['count']
            else:
                nhej_count += g['count']

    wt_count     = type_counts.get('WT', 0)
    edited_count = aligned - wt_count

    # ── Colour palette (professional light theme) ────────────────────────────
    PALETTE = {
        'WT':    '#16a34a',  # professional green
        'DEL':   '#dc2626',  # red
        'INS':   '#2563eb',  # blue
        'SNP':   '#d97706',  # amber
        'KI':    '#7c3aed',  # violet
        'PE':    '#0891b2',  # teal
        'BE':    '#d97706',  # amber
        'bg':    '#ffffff',  # white page background
        'panel': '#f1f5f9',  # very light grey panel
        'text':  '#1e293b',  # near-black body text
        'muted': '#64748b',  # slate grey secondary text
        'accent':'#1e3a5f',  # deep navy (titles / accent)
        'border':'#e2e8f0',  # light border / spine colour
    }

    # ── Base colours (readable on white) ─────────────────────────────────────
    BASE_COLOURS = {'A': '#16a34a', 'T': '#dc2626', 'C': '#1d4ed8', 'G': '#b45309', 'N': '#9ca3af'}

    def base_colour(ch):
        return BASE_COLOURS.get(ch.upper(), '#888888')

    # ── Helpers ──────────────────────────────────────────────────────────────
    def style_panel(ax, colour=None):
        """Style axis as a clean light panel."""
        c = colour or PALETTE['panel']
        ax.set_facecolor(c)
        for spine in ax.spines.values():
            spine.set_color(PALETTE['border'])

    def render_alignment_block(ax, ref_aln, read_aln, title, pct, count,
                               classification, microhomology=None):
        """
        Draw a compact, colour-coded alignment block inside *ax*.
        White-page theme: pastel highlight backgrounds, dark base letters.
        """
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        style_panel(ax)

        c_type    = classification.get('type', '?')
        delta     = classification.get('delta', 0)
        frame     = classification.get('frame', '-')
        delta_str = f"{delta:+d}" if delta != 0 else "±0"
        type_col  = PALETTE.get(c_type, '#94a3b8')

        mh_str = ''
        if microhomology and microhomology.get('length', 0) >= 2:
            mh_str = f"  [{microhomology['pathway']} MH:{microhomology['length']}bp]"

        # ── Left accent bar (type colour) ─────────────────────────────────────
        accent_bar = FancyBboxPatch((0, 0), 0.006, 1,
                                    boxstyle="square,pad=0",
                                    facecolor=type_col, edgecolor='none',
                                    transform=ax.transAxes, clip_on=True, zorder=2)
        ax.add_patch(accent_bar)

        # ── Header ───────────────────────────────────────────────────────────
        header = (f"{title}   {c_type} {delta_str}   {frame}{mh_str}"
                  f"   {count:,} reads  ({pct:.1f}%)")
        ax.text(0.012, 0.88, header, transform=ax.transAxes,
                fontsize=8.5, color=type_col, fontweight='bold',
                fontfamily='monospace', va='center')

        # Thin separator below header
        ax.plot([0.009, 0.999], [0.78, 0.78], transform=ax.transAxes,
                color=PALETTE['border'], linewidth=0.6, clip_on=False)

        # ── Alignment rows (NO CHUNKING — display full sequence once) ─────────
        ROW_SEP = 0.095         # vertical gap between Ref and Rd (8pt text)
        MARKER_SEP = 0.055      # gap from Rd to marker row

        x_offset = 0.025        # space for row label
        available_width = 1.0 - x_offset - 0.01  # leave margin on right
        cw = available_width / len(ref_aln)      # char width scales with seq length
        HALF_H = 0.033          # half-height of highlight rect

        y_ref = 0.68            # Ref row y position
        y_read = y_ref - ROW_SEP  # Read row y position
        y_marker = y_read - MARKER_SEP  # Marker row y position

        # ── Ref row ──
        ax.text(0.008, y_ref, 'Ref:', transform=ax.transAxes,
                fontsize=8, color=PALETTE['muted'],
                fontfamily='monospace', va='center', ha='right',
                clip_on=False)

        for col_idx, ch in enumerate(ref_aln):
            x    = x_offset + col_idx * cw
            r_ch = ref_aln[col_idx]
            q_ch = read_aln[col_idx]

            if r_ch == '-':
                # Insertion in read
                bg = '#dbeafe'; edge_col = '#1d4ed8'
                fc = base_colour(q_ch) if q_ch != '-' else '#9ca3af'
            elif q_ch == '-':
                # Deletion in read
                bg = None; fc = base_colour(ch)
            elif r_ch == q_ch:
                # Match
                bg = None; fc = base_colour(ch)
            else:
                # Mismatch
                bg = '#fef9c3'; edge_col = '#d97706'; fc = '#78350f'

            if bg:
                rect = FancyBboxPatch((x, y_ref - HALF_H), cw * 0.9, HALF_H * 2,
                                     boxstyle="square,pad=0",
                                     facecolor=bg, edgecolor=edge_col, linewidth=0.5,
                                     transform=ax.transAxes, clip_on=True)
                ax.add_patch(rect)

            ax.text(x + cw * 0.5, y_ref, ch,
                    transform=ax.transAxes,
                    fontsize=8, color=fc,
                    fontfamily='monospace', va='center', ha='center')

        # ── Read row ──
        ax.text(0.008, y_read, 'Rd:', transform=ax.transAxes,
                fontsize=8, color=PALETTE['muted'],
                fontfamily='monospace', va='center', ha='right',
                clip_on=False)

        for col_idx, ch in enumerate(read_aln):
            x    = x_offset + col_idx * cw
            r_ch = ref_aln[col_idx]
            q_ch = read_aln[col_idx]

            if r_ch == '-':
                # Insertion in read
                bg = '#dbeafe'; edge_col = '#1d4ed8'
                fc = base_colour(q_ch) if q_ch != '-' else '#9ca3af'
            elif q_ch == '-':
                # Deletion in read
                bg = '#fee2e2'; edge_col = '#dc2626'; fc = '#7f1d1d'
            elif r_ch == q_ch:
                # Match
                bg = None; fc = base_colour(ch)
            else:
                # Mismatch
                bg = '#fef9c3'; edge_col = '#d97706'; fc = '#78350f'

            if bg:
                rect = FancyBboxPatch((x, y_read - HALF_H), cw * 0.9, HALF_H * 2,
                                     boxstyle="square,pad=0",
                                     facecolor=bg, edgecolor=edge_col, linewidth=0.5,
                                     transform=ax.transAxes, clip_on=True)
                ax.add_patch(rect)

            ax.text(x + cw * 0.5, y_read, ch,
                    transform=ax.transAxes,
                    fontsize=8, color=fc,
                    fontfamily='monospace', va='center', ha='center')

        # ── Marker row (shows differences) ──
        for col_idx in range(len(ref_aln)):
            x    = x_offset + col_idx * cw
            r_ch = ref_aln[col_idx]
            q_ch = read_aln[col_idx]

            if r_ch != q_ch or r_ch == '-' or q_ch == '-':
                ax.text(x + cw * 0.5, y_marker, '^',
                        transform=ax.transAxes,
                        fontsize=6, color=PALETTE['muted'],
                        fontfamily='monospace', va='center', ha='center')

    # ════════════════════════════════════════════════════════════════════════
    # PDF PAGES
    # ════════════════════════════════════════════════════════════════════════
    now_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    with PdfPages(output_path) as pdf:

        # ════════════════════════════════════════════════════════════════
        # PAGE 1 — SUMMARY DASHBOARD
        # ════════════════════════════════════════════════════════════════
        fig = plt.figure(figsize=(11.69, 8.27), facecolor=PALETTE['bg'])  # A4 landscape
        gs  = gridspec.GridSpec(2, 4, figure=fig,
                                left=0.06, right=0.97,
                                top=0.82,  bottom=0.10,
                                hspace=0.50, wspace=0.45)

        # ── Decorative header band ────────────────────────────────────────
        header_ax = fig.add_axes([0, 0.964, 1, 0.036])
        header_ax.set_facecolor(PALETTE['accent'])
        header_ax.axis('off')
        fig.add_artist(plt.Line2D([0, 1], [0.962, 0.962],
                                  transform=fig.transFigure,
                                  color=PALETTE['border'], linewidth=0.6))

        # ── Page header ──────────────────────────────────────────────────
        fig.text(0.015, 0.978, 'CleanFinder', fontsize=14, color='#ffffff',
                 fontweight='bold', fontfamily='monospace', va='center')
        fig.text(0.20, 0.978,
                 f'Sample: {sample_name}   ·   Mode: {mode.upper()}   ·   {now_str}',
                 fontsize=8.5, color='#cbd5e1', va='center')
        fig.text(0.97, 0.978,
                 'IUF Leibniz Research Institute, Düsseldorf',
                 fontsize=7.5, color='#94a3b8', va='center', ha='right')
        # ── Summary statistics table ─────────────────────────────────────────
        in_f = frame_counts['In-Frame']
        out_f = frame_counts['Out-Frame']
        edit_pct = 100 * edited_count / max(aligned, 1)
        in_frame_pct = 100 * in_f / max(in_f + out_f, 1) if (in_f + out_f) > 0 else 0

        # Draw clean summary table
        summary_y = 0.935
        col1_x, col2_x, col3_x = 0.06, 0.35, 0.64
        line_h = 0.0165

        summary_data = [
            ('Sequencing', '', ''),
            ('Total Reads', f'{total:,}', f'({unaligned:,} unaligned)'),
            ('Aligned', f'{aligned:,}', f'({100*aligned/max(total,1):.1f}%)'),
            ('', '', ''),
            ('Editing', '', ''),
            ('Wildtype', f'{wt_count:,}', f'({100*wt_count/max(aligned,1):.1f}%)'),
            ('Edited', f'{edited_count:,}', f'({edit_pct:.1f}%)'),
            ('', '', ''),
            ('Frame Status (edits only)', '', ''),
            ('In-Frame', f'{in_f:,}', f'({in_frame_pct:.1f}%)'),
            ('Out-of-Frame', f'{out_f:,}', f'({100-in_frame_pct:.1f}%)'),
        ]

        for label, val1, val2 in summary_data:
            if label == '':
                summary_y -= line_h * 0.5
                continue

            # Category headers
            if val1 == '' and val2 == '':
                fig.text(col1_x, summary_y, label, fontsize=8, color=PALETTE['accent'],
                        fontweight='bold', va='top')
            else:
                fig.text(col1_x, summary_y, label, fontsize=7.5,
                        color=PALETTE['muted'], va='top')
                fig.text(col2_x, summary_y, val1, fontsize=7.5,
                        color=PALETTE['text'], va='top', fontfamily='monospace', fontweight='bold')
                fig.text(col3_x, summary_y, val2, fontsize=7.5,
                        color=PALETTE['muted'], va='top', fontfamily='monospace')

            summary_y -= line_h

        # ── Allele list table ─────────────────────────────────────────────────
        allele_y = summary_y - 0.025  # Start after summary with a small gap
        col_rank_x = 0.06
        col_name_x = 0.12
        col_count_x = 0.40
        col_pct_x = 0.55
        col_frame_x = 0.70
        line_h = 0.016

        # Header
        fig.text(col_rank_x, allele_y, '#', fontsize=8, color=PALETTE['accent'],
                fontweight='bold', va='top')
        fig.text(col_name_x, allele_y, 'Allele', fontsize=8, color=PALETTE['accent'],
                fontweight='bold', va='top')
        fig.text(col_count_x, allele_y, 'Reads', fontsize=8, color=PALETTE['accent'],
                fontweight='bold', va='top', ha='right')
        fig.text(col_pct_x, allele_y, 'Percentage', fontsize=8, color=PALETTE['accent'],
                fontweight='bold', va='top', ha='right')
        fig.text(col_frame_x, allele_y, 'Frame', fontsize=8, color=PALETTE['accent'],
                fontweight='bold', va='top')
        allele_y -= line_h * 1.2

        # Draw alleles (all if show_all, else top 15)
        max_alleles = len(sorted_genos) if show_all else 15
        for rank, (label, geno) in enumerate(sorted_genos[:max_alleles], 1):
            allele_type = geno['classification']['type']
            delta = geno['classification'].get('delta', 0)
            count = geno['count']
            frame = geno['classification'].get('frame', '—')
            pct = 100 * count / max(aligned, 1)

            # Create label
            if allele_type == 'WT':
                allele_label = 'WT ±0'
            elif allele_type == 'DEL':
                allele_label = f'DEL {delta}'
            elif allele_type == 'INS':
                allele_label = f'INS {delta:+d}'
            else:
                allele_label = f'{allele_type}'

            # Color
            if allele_type == 'WT':
                col = PALETTE['WT']
            elif allele_type == 'DEL':
                intensity = min(abs(delta) / 30.0, 1.0)
                r = int(220 - intensity * 50)
                g = int(40 + intensity * 20)
                b = int(40 + intensity * 20)
                col = f'#{r:02x}{g:02x}{b:02x}'
            elif allele_type == 'INS':
                col = PALETTE['INS']
            else:
                col = PALETTE.get(allele_type, '#94a3b8')

            fig.text(col_rank_x, allele_y, str(rank), fontsize=7.5,
                    color=PALETTE['muted'], va='top')
            fig.text(col_name_x, allele_y, allele_label, fontsize=7.5,
                    color=col, fontweight='bold', va='top')
            fig.text(col_count_x, allele_y, f'{count:,}', fontsize=7.5,
                    color=PALETTE['text'], va='top', ha='right', fontfamily='monospace')
            fig.text(col_pct_x, allele_y, f'{pct:.1f}%', fontsize=7.5,
                    color=PALETTE['text'], va='top', ha='right', fontfamily='monospace')
            fig.text(col_frame_x, allele_y, frame, fontsize=7, color=PALETTE['muted'],
                    va='top')

            allele_y -= line_h

        # ── Bottom panel: indel size spectrum ────────────────────────────
        ax5 = fig.add_subplot(gs[1, :])
        style_panel(ax5)
        ax5.yaxis.grid(True, linestyle=':', linewidth=0.5, color=PALETTE['border'], zorder=0)
        ax5.set_axisbelow(True)

        # Collect delta values for DEL and INS
        delta_counts = defaultdict(int)
        for _, g in sorted_genos:
            c = g['classification']
            if c['type'] in ('DEL', 'INS') and c['delta'] != 0:
                delta_counts[c['delta']] += g['count']

        if delta_counts:
            deltas = sorted(delta_counts.keys())
            counts = [delta_counts[d] for d in deltas]
            bar_cols = [PALETTE['DEL'] if d < 0 else PALETTE['INS'] for d in deltas]
            x_pos = range(len(deltas))
            ax5.bar(x_pos, counts, color=bar_cols, alpha=0.85, width=0.7)
            ax5.set_xticks(list(x_pos))
            ax5.set_xticklabels([str(d) for d in deltas],
                                fontsize=6.5, color=PALETTE['muted'], rotation=45)
            ax5.tick_params(axis='y', colors=PALETTE['muted'], labelsize=7)
            ax5.set_xlabel('Indel size (bp)',  color=PALETTE['muted'], fontsize=8)
            ax5.set_ylabel('Read count',       color=PALETTE['muted'], fontsize=8)

            legend_p = [mpatches.Patch(color=PALETTE['DEL'], label='Deletion'),
                        mpatches.Patch(color=PALETTE['INS'], label='Insertion')]
            ax5.legend(handles=legend_p, fontsize=7.5,
                       facecolor=PALETTE['panel'], edgecolor=PALETTE['border'],
                       labelcolor=PALETTE['text'])
        else:
            ax5.text(0.5, 0.5, 'No indels detected', ha='center', va='center',
                     color=PALETTE['muted'], fontsize=11, transform=ax5.transAxes)

        ax5.set_title('Indel Size Spectrum', color=PALETTE['text'], fontsize=9, pad=6)
        for spine in ax5.spines.values():
            spine.set_color(PALETTE['border'])

        # Footer
        fig.text(0.97, 0.015, 'CleanFinder  ·  IUF Leibniz Research Institute, Düsseldorf',
                 ha='right', fontsize=7, color=PALETTE['border'])

        pdf.savefig(fig, facecolor=fig.get_facecolor())
        plt.close(fig)

        # ════════════════════════════════════════════════════════════════
        # PAGES 2+ — TOP ALIGNMENTS  (5 per page, precisely positioned)
        # ════════════════════════════════════════════════════════════════
        ALNS_PER_PAGE = 5
        BLOCK_H  = 0.156   # fraction of figure height per alignment block
        BLOCK_GAP = 0.007  # gap between blocks
        BLK_LEFT, BLK_W = 0.025, 0.968  # horizontal extent — wider, less wasted margin

        top_genos = [(sig, g) for sig, g in sorted_genos
                     if g['alignment'] and g['alignment']['ok']]
        if not show_all:
            top_genos = top_genos[:15]

        for page_start in range(0, len(top_genos), ALNS_PER_PAGE):
            batch = top_genos[page_start:page_start + ALNS_PER_PAGE]

            fig2 = plt.figure(figsize=(11.69, 8.27), facecolor=PALETTE['bg'])
            # ── Decorative header band ────────────────────────────────────
            hax2 = fig2.add_axes([0, 0.964, 1, 0.036])
            hax2.set_facecolor(PALETTE['accent'])
            hax2.axis('off')
            fig2.add_artist(plt.Line2D([0, 1], [0.962, 0.962],
                                       transform=fig2.transFigure,
                                       color=PALETTE['border'], linewidth=0.6))
            fig2.text(0.015, 0.978, 'CleanFinder — Top Genotype Alignments',
                      fontsize=13, color='#ffffff',
                      fontweight='bold', fontfamily='monospace', va='center')
            fig2.text(0.60, 0.978, f'Sample: {sample_name}   ·   Mode: {mode.upper()}',
                      fontsize=8.5, color='#cbd5e1', va='center')

            # ── Place each block with exact coordinates ───────────────────
            y_cursor = 0.953  # start just below the header band

            for slot, (sig, g) in enumerate(batch):
                rank  = page_start + slot + 1
                count = g['count']
                pct   = 100 * count / max(aligned, 1)
                aln   = g['alignment']
                cls   = g['classification']
                mh    = g.get('microhomology')

                y_bottom = y_cursor - BLOCK_H
                ax_aln   = fig2.add_axes([BLK_LEFT, y_bottom, BLK_W, BLOCK_H])
                render_alignment_block(ax_aln, aln['alnRef'], aln['alnRead'],
                                       f'#{rank}', pct, count, cls, mh)
                y_cursor = y_bottom - BLOCK_GAP

            # ── Color legend (bottom of page) ─────────────────────────────
            leg_y = max(y_cursor - 0.002, 0.018)
            legend_items = [
                ('#dbeafe', '#1d4ed8', 'Insertion'),
                ('#fee2e2', '#7f1d1d', 'Deletion'),
                ('#fef9c3', '#78350f', 'Mismatch / SNP'),
            ]
            lx = BLK_LEFT
            for bg_c, txt_c, label in legend_items:
                swatch = FancyBboxPatch((lx, leg_y), 0.018, 0.012,
                                        boxstyle="square,pad=0",
                                        facecolor=bg_c, edgecolor=PALETTE['border'],
                                        linewidth=0.5,
                                        transform=fig2.transFigure, clip_on=False)
                fig2.add_artist(swatch)
                fig2.text(lx + 0.022, leg_y + 0.006, label,
                          fontsize=6.5, color=PALETTE['muted'], va='center',
                          transform=fig2.transFigure)
                lx += 0.11

            fig2.text(0.97, 0.015,
                      'CleanFinder  ·  IUF Leibniz Research Institute, Düsseldorf',
                      ha='right', fontsize=7, color=PALETTE['border'])

            pdf.savefig(fig2, facecolor=fig2.get_facecolor())
            plt.close(fig2)

    print(f"\n📄 PDF report saved to: {output_path}")


if __name__ == '__main__':
    main()
