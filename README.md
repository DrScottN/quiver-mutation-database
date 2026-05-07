# Mathematicians vs Machines: 
## Baseline Algorithms for Mutation Equivalence and Acyclicity
This is the code repository for our paper, of the same name.
It provides two main things (and several minor ones):
1. A collection of baseline algorithms for quiver mutation-equivalence and mutation-acyclicity (available through `prediction.py`);
2. the code to generate additional mutation-equivalence and mutation-acyclicity datasets with many associated statistics and labels (available through `generation.py`).

Both of these tasks rely heavily on our self-contained implementation of quiver mutation, mutation classes, and associated invariants (in `quiver.py`).


## Requirements

`pip install -r requirements.txt`

This repository uses:

`Python 3.12+`
`numpy`

Optionally, for `test_croissant.py`:
`mlcroissant`


## Databset

Our Dataset is also available on [Huggingface](https://huggingface.co/datasets/labeledquivers/LabeledQuiverMutationAcyclicity).

The dataset can be generated from scratch via:
`python make4VertexMA.py`

## Baselines

We provide algorithms based on known mutation invariants for mutation equivalence and detecting mutation acyclicity. We evaluate these algorithms on many previously considered families of quivers, and find they do extremely well. 

Run `python exam.py` to repeat these experiments. See the code therein for example useage.

## Additional Contents

### Mutation Acyclic Results
This folder contains print-outs from calling `python generation.py` with a range of settings.
These were used primarily in estimating runtime, debugging, and profiling.

### Model Notebooks
These contain the google-colab python notebooks used for our trained models on the dataset. Based heavily on examples from a workshop by [Kyu-Hwan Lee](https://khlee-math.github.io/).

### CSV scripts
Both `deduplication.py` and `difference.py` are cleanup and debugging scripts for resulting csv files. `deduplication.py` reads a csv indexed by an exchange matrix, with some ternary fields and potential duplicates, and combines duplicated entries while preserving ternary information.
`difference.py` takes the set difference of two csv files and writes it to a new csv.

### polynomials
We need integer polynomials to compute the Alexander polynomial. We use the polynomial class from `polynomial.py`. This has its own set of tests in `test_polys.py`.

### unknown
We implement ternary logic using the unknown class, which overrides the logical operators `& |`. The `Unknown()` object functions like a `NULL` entry in a SQL table. 

## Testing
We provide a suite of tests for our implementation of quivers and their mutation-invariants, available in `test.py`.
To run our tests, use
`python test.py`
from the directory.


## Contribute
We welcome contributions to this code! You may reach us via creating an issue here on github.