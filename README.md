# Phylogenetic Tree Construction and Evaluation

This project provides a streamlined pipeline for building and evaluating phylogenetic trees from aligned sequence data. It supports tree construction via multiple algorithms and compares their output using established metrics. The system is batch-ready and can process multiple FASTA files in one run, skipping invalid or misaligned inputs.

---

## Project Structure

```
.
├── input/                         # Input folder with aligned FASTA files
│   ├── file1.fasta
│   ├── file2.fasta
│   └── ...
├── output_images/                # Tree visualizations organized per FASTA file
│   ├── file1/
│   │   ├── UPGMA.png
│   │   └── ...
│   └── file2/
│       └── ...
├── output/                       # Text summaries for each file
│   ├── file1_output.txt
│   └── ...
├── phylogenetic_tree.py          # Main tree construction and evaluation script
├── requirements.txt              # Dependency list
└── README.md                     # Project documentation
```

---

## Features

- **Batch processing**: Automatically scans the `input/` directory for aligned `.fasta` or `.fa` files
- **Tree construction using**:
  - UPGMA
  - Neighbor Joining
  - Maximum Parsimony
  - Maximum Likelihood *(simulated)*
  - Bayesian Inference *(simulated)*
- **Performance metrics**:
  - Execution Time
  - Memory Usage
  - Robinson-Foulds Distance
  - Symmetric Branch Difference
- **Output handling**:
  - Creates `output_images/<fasta_file>/` and saves tree images
  - Creates `output/<fasta_file>_output.txt` with performance and comparison results
- **Graceful error handling**: Skips files that fail alignment or tree building

---

## Installation

Make sure Python 3.8+ is installed. Then install the required packages:

```bash
pip install -r requirements.txt
```

---

## Usage

Place all your aligned FASTA files in the `input/` directory:

```
input/
├── aligned1.fasta
├── aligned2.fasta
└── ...
```

Then run the script:

```bash
python phylogenetic_tree.py
```

Each valid input will generate:

- A subfolder in `output_images/` with PNG visualizations for all tree methods
- A metrics report in `output/` summarizing:
  - Time
  - Memory
  - RF Distance
  - KC Distance

---

## Example Console Output

```
Processing file: aligned1
Evaluating UPGMA...
UPGMA: Time = 0.001234s | Mem = 32.00 KB | RF = 0 | KC = 0
Saved output_images/aligned1/UPGMA.png

Evaluating Neighbor Joining...
Neighbor Joining: Time = 0.001562s | Mem = 28.00 KB | RF = 2 | KC = 3
...

Saved output/aligned1_output.txt
```

---

## Notes

- **Input Files**: Must be aligned FASTA format. You can use tools like MAFFT or MUSCLE to prepare inputs.
- **Simulated Models**: Maximum Likelihood and Bayesian Inference are placeholders and can be replaced with real implementations.
- **Evaluation Reference**: UPGMA is used as the baseline tree for all comparison metrics.

---

## License

This project is licensed under the MIT License.
