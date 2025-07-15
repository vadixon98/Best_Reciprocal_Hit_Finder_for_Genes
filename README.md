# Best Reciprocal Hit (BRH) Finder for Gene Homology

A Python-based toolset to identify **Best Reciprocal Hits (BRHs)** between genes from two species, using protein sequence similarity scored by the **BLOSUM62** substitution matrix and dynamic programming alignment.

---

## Requirements

* **Python 3.x**
* No external dependencies (standard library only).

---

## Installation

1. Clone or download the repository:

   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```
2. (Optional) Create and activate a virtual environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # Windows: venv\Scripts\activate
   ```

---

## Scripts Overview

### 1. `Best Reciprocal Hit Finder for Genes.py`

Implements the core BRH workflow:

* `alignScore(seq1, seq2, gap, matrix)`: Recursive Needleman–Wunsch alignment score.
* `memoAlignScore(seq1, seq2, gap, matrix, memo)`: Memoized version for speed.
* `allScores(list1, list2)`: Compute pairwise alignment scores for two gene lists.
* `closestMatch(gene, scoresDict)`: Identify the highest-scoring partner in the other species.
* `printBRH(gene, scoresDict)`: Print BRH pairs when reciprocal.
* `runBRHSample()`: Execute sample analysis on small predefined datasets.
* `runBRH()`: Run full analysis on human vs. chicken datasets.

**Usage example**:

```python
from Best_Reciprocal_Hit_Finder_for_Genes import runBRHSample
runBRHSample()
```

### 2. `blosum62.py`

Provides the **BLOSUM62** substitution matrix as a Python dictionary for scoring amino acid substitutions.

### 3. `humanChickenProteins.py`

Contains gene metadata and sequence lists:

* `humanGeneList` / `chickenGeneList`: Full species datasets.
* `sampleHumanGeneList` / `sampleChickenGeneList`: Smaller test sets.
* `geneD`: Dictionary with gene identifiers and corresponding protein sequences.

---

## Usage Tips

* **Run full BRH analysis**:

  ```python
  from Best_Reciprocal_Hit_Finder_for_Genes import runBRH
  runBRH()
  ```
* **Adjust gap penalties** by passing a different `gap` value to `alignScore` or `memoAlignScore`.
* **Optimize performance** by using the memoized alignment function for larger datasets.
* **Customize datasets** by editing `humanChickenProteins.py` or supplying your own lists of sequences.

---

## Contributing

Contributions and issues are welcome! Please fork the repository, make your changes, and submit a pull request on GitHub.

---

## License

This project is released under the MIT License. Include a copy of the license if distributing.
