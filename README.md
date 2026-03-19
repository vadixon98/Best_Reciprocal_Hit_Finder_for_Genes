<div align="center">

# Best Reciprocal Hit (BRH) Finder for Genes

*A Python tool for identifying Best Reciprocal Hits between genes from different species using protein sequence alignment and the BLOSUM62 substitution matrix.*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Comparative%20Genomics-orange)](https://github.com)
[![Algorithm](https://img.shields.io/badge/Algorithm-Dynamic%20Programming-purple)](https://github.com)

---

</div>

## Features

- **Optimal Alignment Scoring** — Needleman–Wunsch algorithm with BLOSUM62 substitution matrix
- **Performance Optimized** — Memoized dynamic programming for efficient computation on larger datasets
- **Reciprocal Hit Detection** — Identifies bidirectional best matches (orthologs)
- **Pairwise Analysis** — Compares all genes between two species
- **Sample and Full Datasets** — Includes sample data for testing and full human vs. chicken datasets
- **Documented Code** — Functions include docstrings and inline comments

---

## Quick Start

### Prerequisites

- Python 3.8+ (Python 3.10+ recommended)
- Standard library only; no external dependencies

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-username/Best_Reciprocal_Hit_Finder_for_Genes.git
   cd Best_Reciprocal_Hit_Finder_for_Genes
   ```

2. **(Optional) Create a virtual environment:**
   ```bash
   # Windows
   python -m venv venv
   venv\Scripts\activate
   
   # macOS/Linux
   python3 -m venv venv
   source venv/bin/activate
   ```

3. No additional packages are required.

---

## What is a Best Reciprocal Hit?

A **Best Reciprocal Hit (BRH)** occurs when:

> Gene A from species 1 matches most closely with Gene B from species 2  
> **AND**  
> Gene B from species 2 matches most closely with Gene A from species 1

This bidirectional relationship is a strong indicator of **orthology** — genes that descended from a common ancestral gene and typically retain similar functions across species.

```
Species 1          Species 2
  Gene A   ════════>   Gene B  (A's best match)
  Gene A   <════════   Gene B  (B's best match)
              BRH
```

---

## Scripts Overview

### Main Script: `best_reciprocal_hit_finder.py`

The core BRH analysis toolkit:

| Function | Description | Use Case |
|----------|-------------|----------|
| `alignScore()` | Recursive Needleman–Wunsch alignment | Basic alignment scoring |
| `memoAlignScore()` | Memoized alignment | Large sequence datasets |
| `allScores()` | Pairwise scores for all gene pairs | Build score matrix |
| `closestMatch()` | Find highest-scoring partner | Identify best match |
| `printBRH()` | Print reciprocal hit pairs | Display BRH results |
| `runBRHSample()` | Run on sample datasets | Quick testing |
| `runBRH()` | Full human vs. chicken analysis | Complete analysis |

### Supporting Modules

#### `blosum62.py`
Provides the **BLOSUM62** substitution matrix as a Python dictionary for scoring amino acid substitutions during alignment.

#### `humanChickenProteins.py`
Contains gene metadata and sequence data:
- `humanGeneList` / `chickenGeneList` → Full species datasets
- `sampleHumanGeneList` / `sampleChickenGeneList` → Smaller test sets
- `geneD` → Dictionary mapping gene IDs to protein sequences

---

## Usage Examples

### Sample Data

```python
from best_reciprocal_hit_finder import runBRHSample

# Run analysis on sample datasets
runBRHSample()
```

**Output:**
```
human --- chicken
chr1 100000 GeneA --- chr2 50000 GeneX
chr1 200000 GeneB --- chr5 75000 GeneY
...

chicken --- human
chr2 50000 GeneX --- chr1 100000 GeneA
...
```

### Full Dataset Analysis

```python
from best_reciprocal_hit_finder import runBRH

# ⚠️ This may take considerable time for large datasets
runBRH()
```

### Custom Alignment

```python
from best_reciprocal_hit_finder import memoAlignScore
from blosum62 import blosum62

# Compute alignment score between two sequences
memo = {}
score = memoAlignScore(
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
    "MKTIIALSYIFCLVFAADPERKYLVEARARLERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL",
    -9,
    blosum62,
    memo
)
print(f"Alignment score: {score}")
```

---

## Customization

### Gap Penalties

```python
# Default gap penalty is -9
# You can modify it for different sensitivity
score = memoAlignScore(seq1, seq2, gap=-10, substitutionMatrix=blosum62, memo={})
```

### Performance

- Use `memoAlignScore()` for large datasets
- The shared `memo` dictionary in `allScores()` reuses cached subproblem solutions
- For very large analyses, consider batch processing or parallelization

### Custom Datasets

All genes must exist in `geneD` (from `humanChickenProteins.py`) with their protein sequences. Use a subset of the built-in genes, or add your own data to `humanChickenProteins.py` first.

```python
# Use a subset of genes from the built-in dataset
# All genes must exist in geneD
from best_reciprocal_hit_finder import allScores, printBRH

# Subsets of the sample gene lists (h4, h6, c8, etc. are in geneD)
myHumanGenes = ['h4', 'h6', 'h9']
myChickenGenes = ['c8', 'c17', 'c19']

allScoresD = allScores(myHumanGenes, myChickenGenes)

for gene in myHumanGenes:
    printBRH(gene, allScoresD)
```

---

## Algorithm

### Alignment

The tool uses the **Needleman-Wunsch** global alignment algorithm with three operations at each position:

1. **Match/Mismatch** → Align two amino acids (scored by BLOSUM62)
2. **Gap in Sequence 1** → Skip character in S1
3. **Gap in Sequence 2** → Skip character in S2

The algorithm recursively finds the optimal alignment by choosing the maximum score among these three options, with memoization to cache intermediate results for efficiency.

---

## Contributing

Contributions are welcome. To contribute:

1. Fork the repository
2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/Best_Reciprocal_Hit_Finder_for_Genes.git
   ```
3. Create a branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```
4. Make your changes and test them
5. Commit with clear messages
6. Push to your fork and open a Pull Request

---

## License

This project is released under the MIT License. See the license file for details.
