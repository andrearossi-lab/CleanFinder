

 ____ _                  _____ _           _           
  / ___| | ___  __ _ _ __ |  ___(_)_ __   __| | ___ _ __ 
 | |   | |/ _ \/ _` | '_ \| |_  | | '_ \ / _` |/ _ \ '__|
 | |___| |  __/ (_| | | | |  _| | | | | | (_| |  __/ |   
  \____|_|\___|\__,_|_| |_|_|   |_|_| |_|\__,_|\___|_|  

Algorithm

CleanFinder uses a two-stage, anchor-constrained glocal alignment strategy to classify CRISPR editing outcomes.

Stage 1: Anchor-Based Read Extraction

Each read is screened for two flanking anchor sequences (default: 20 bp) derived from the reference edges around the gRNA cut site. Anchors are matched using a fuzzy search with configurable mismatch tolerance, supporting both forward and reverse-complement orientations.

```
Reference:   [── Left Anchor (20bp) ──][── Analysis Window ──][── Right Anchor (20bp) ──]
                                        ↑ gRNA cut site

Read:        ...NNNN[── L Anchor ──][── Extracted Middle ──][── R Anchor ──]NNNN...
```

Reads that fail to match both anchors are discarded as off-target. An optional **K-mer pre-filter** (10–13 bp seed from the left anchor) provides early rejection of unrelated reads for high-throughput performance.

Stage 2: Glocal Fitting Alignment

The extracted middle region is aligned against the reference window using a **glocal (fitting) alignment** with **affine gap penalties**. This forces the entire reference to align end-to-end against a substring of the read, naturally accommodating insertions and deletions at the cut site.

The dynamic programming engine uses three state matrices (M, Ix, Iy) stored as flat `Int32Array` buffers for memory-efficient computation:

```javascript
// Core alignment: 3-state affine gap DP
// M  = Match/Mismatch matrix
// Ix = Gap in Read (Deletion relative to Reference)
// Iy = Gap in Reference (Insertion relative to Reference)

function fittingAlign(ref, read, match = 2, mismatch = -6, gapOpen = -10, gapExt = -2) {
    const m = ref.length, n = read.length;
    const width = n + 1;
    const size = (m + 1) * width;
    const M  = new Int32Array(size);
    const Ix = new Int32Array(size);
    const Iy = new Int32Array(size);

    // Initialization: free start in read (fitting), penalized gaps in ref (global)
    M[0] = 0;
    for (let j = 1; j <= n; j++) M[j] = 0;           // Free read start
    for (let i = 1; i <= m; i++) Ix[i * width] = gapOpen + (i - 1) * gapExt;

    // Fill DP
    for (let i = 1; i <= m; i++) {
        for (let j = 1; j <= n; j++) {
            const cur = i * width + j;
            const score = (ref[i-1] === read[j-1]) ? match : mismatch;

            Iy[cur] = Math.max(M[cur-1] + gapOpen, Iy[cur-1] + gapExt);
            Ix[cur] = Math.max(M[(i-1)*width+j] + gapOpen, Ix[(i-1)*width+j] + gapExt);

            const diag = (i-1) * width + (j-1);
            M[cur] = Math.max(M[diag], Ix[diag], Iy[diag]) + score;
        }
    }
    // Traceback from last row (ref fully consumed) → free read end
    // ...
}
```

Mode-specific scoring presets:

| Mode | Match | Mismatch | Gap Open | Gap Extend | Rationale |
|------|-------|----------|----------|------------|-----------|
| Cas9 | +2 | −6 | −3 | −1 | Balanced for short indels |
| Cas12 | +2 | −6 | −5 | −0.5 | Accommodates staggered cuts |
| Prime Editing | +2 | −3 | −1 | −0.3 | Permissive for long insertions |
| Base Editing | +2 | −4 | −20 | −2 | Prohibits indels, favors substitutions |

Genotype Classification

After alignment, each read is classified based on the length difference (Δ) between extracted middle and reference:

- **Δ = 0, identical** → Wild-Type (WT)
- **Δ = 0, mismatches** → SNP or Base Edit
- **Δ < 0** → Deletion (In-Frame if Δ mod 3 = 0)
- **Δ > 0** → Insertion (In-Frame if Δ mod 3 = 0)
- **Template match** → Knock-in (HDR) or Prime Edit

Deletions are further analyzed for **microhomology** at junction sites to classify repair pathway (MMEJ ≥ 2 bp vs. NHEJ).

Project Structure

```
CleanFinder/
├── index_v5.html           # Web application (main entry point)
├── cleanfinder.py          # Python CLI for batch/pipeline use
├── js/
│   └── app_v5.js           # Core engine: alignment, analysis, UI
├── css/
│   └── style.css           # Application styles
├── turbomode.html          # Turbo Mode (high-throughput screens)
├── allelic_dropout.html    # Allelic Dropout detection (beta)
├── genome_viewer.html      # Interactive Genome Viewer
├── README.md
└── LICENSE
```

Citation

> **CleanFinder**: A cross-platform genotyping suite for CRISPR editing analysis.  
> Rossi A. et al., 2026.  
> *Manuscript in preparation.*

License

MIT License — see [LICENSE](LICENSE) for details.

Copyright (c) 2026 Andrea Rossi

Developed by **GEMD** and **EACR** labs, IUF Leibniz Research Institute, Düsseldorf.  
Funded by DFG (Deutsche Forschungsgemeinschaft).
