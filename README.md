# CleanFinder
CleanFinder: Browser-Native Analysis of CRISPR Editing Outcomes

Alignment Engine

CleanFinder uses a client-side implementation of the Needleman–Wunsch global alignment algorithm with affine gap penalties (scoring: match = +2, mismatch = −1, gap open = −5, gap extend = −1).

This algorithm performs per-read alignment between the extracted reference window (refMiddle) and each sequencing read to accurately characterize indels near CRISPR cut sites.

All computations are executed locally in the browser to ensure data privacy and eliminate the need for server-side processing.
