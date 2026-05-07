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



## Testing
We provide a suite of tests for our implementation of quivers, available in test.py.
To run our tests, use
`python test.py`
from the directory.
The tests themselves may also be of interest, as they provide many examples of how to compute with our code.


## Mutation Acyclic Results
This folder contains print-outs from calling `python generation.py` with a range of settings.
These were used primarily in estimating runtime, debugging, and profiling.