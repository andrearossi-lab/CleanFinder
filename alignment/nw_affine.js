/*
  Needleman–Wunsch global alignment with affine gap penalties
  CleanFinder implementation (Spiessbach et al., 2025)
  Scoring scheme:
      match = +2, mismatch = −1, gap open = −5, gap extend = −1
  This algorithm is used in the FASTQ Analyzer for read classification.
*/

function nwAffineAlign(seq1, seq2, scoring = {match: 2, mismatch: -1, gapOpen: -5, gapExtend: -1}) {
  const { match, mismatch, gapOpen, gapExtend } = scoring;
  const m = seq1.length;
  const n = seq2.length;

  // Matrices for matches (M), insertion in seq1 (Ix), deletion in seq2 (Iy)
  const M  = Array.from({ length: m + 1 }, () => new Float32Array(n + 1));
  const Ix = Array.from({ length: m + 1 }, () => new Float32Array(n + 1));
  const Iy = Array.from({ length: m + 1 }, () => new Float32Array(n + 1));

  const NEG_INF = -1e9;

  // Initialize matrices
  for (let i = 0; i <= m; i++) {
    for (let j = 0; j <= n; j++) {
      M[i][j] = Ix[i][j] = Iy[i][j] = NEG_INF;
    }
  }

  M[0][0] = 0;
  for (let i = 1; i <= m; i++) Iy[i][0] = gapOpen + (i - 1) * gapExtend;
  for (let j = 1; j <= n; j++) Ix[0][j] = gapOpen + (j - 1) * gapExtend;

  // Fill matrices
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const s = (seq1[i - 1] === seq2[j - 1]) ? match : mismatch;

      // Match/mismatch
      M[i][j] = Math.max(
        M[i - 1][j - 1],
        Ix[i - 1][j - 1],
        Iy[i - 1][j - 1]
      ) + s;

      // Gap in seq1 (insertion in seq2)
      Ix[i][j] = Math.max(
        M[i][j - 1] + gapOpen,
        Ix[i][j - 1] + gapExtend
      );

      // Gap in seq2 (deletion in seq1)
      Iy[i][j] = Math.max(
        M[i - 1][j] + gapOpen,
        Iy[i - 1][j] + gapExtend
      );
    }
  }

  // Traceback
  let i = m, j = n;
  let aligned1 = "", aligned2 = "";
  let matrix = "M";
  let score = Math.max(M[i][j], Ix[i][j], Iy[i][j]);

  if (score === Ix[i][j]) matrix = "X";
  if (score === Iy[i][j]) matrix = "Y";

  while (i > 0 || j > 0) {
    if (matrix === "M") {
      if (i > 0 && j > 0) {
        const s = (seq1[i - 1] === seq2[j - 1]) ? match : mismatch;
        if (M[i][j] === M[i - 1][j - 1] + s ||
            M[i][j] === Ix[i - 1][j - 1] + s ||
            M[i][j] === Iy[i - 1][j - 1] + s) {
          aligned1 = seq1[i - 1] + aligned1;
          aligned2 = seq2[j - 1] + aligned2;
          i--; j--;
          continue;
        }
      }
    }
    if (matrix === "X" || (j > 0 && M[i][j] === Ix[i][j])) {
      aligned1 = "-" + aligned1;
      aligned2 = seq2[j - 1] + aligned2;
      j--;
      matrix = (Ix[i][j] + gapExtend === Ix[i][j + 1]) ? "X" : "M";
      continue;
    }
    if (matrix === "Y" || (i > 0 && M[i][j] === Iy[i][j])) {
      aligned1 = seq1[i - 1] + aligned1;
      aligned2 = "-" + aligned2;
      i--;
      matrix = (Iy[i][j] + gapExtend === Iy[i + 1][j]) ? "Y" : "M";
      continue;
    }
    // fallback
    if (i > 0 && j > 0) { i--; j--; }
    else if (i > 0) i--;
    else if (j > 0) j--;
  }

  return { aligned1, aligned2, score };
}
