# quiver-mutation-database
A database of small-weight quivers labeled by their mutation classes, along with heuristic algorithms. Datasets are also available on [Huggingface](https://huggingface.co/datasets/labeledquivers/LabeledQuiverMutationAcyclicity).



## Requirements

`pip install -r requirements.txt`

This repository uses:

`Python 3.12+`
`numpy`

Optionally, for `test_croissant.py`:
`mlcroissant`


## Testing
We provide a suite of tests for our implementation of quivers, available in test.py.
To run our tests, use
`python test.py`
from the directory.
The tests themselves may also be of interest, as they provide many examples of how to compute with our code.


## Mutation Acyclic Results
This folder contains print-outs from calling `python generation.py` with a range of settings.
These were used primarily in estimating runtime, debugging, and profiling.