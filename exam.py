from quiver import *
from prediction import *
from generation import *

def main():
    """ 
    Evaluation of several previous quiver mutation databases using the benchmarks in prediction.
     """
    # ML for AG: A7, X7, and the oriented cycle with weights = 2
    A7 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0), dtype=object).reshape((7,7)))
    X7 = Quiver(np.array((0, 2, -1, 0, 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, -1, 0, 1, -1, 1, -1, 0, 0, -1, 0, 2, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 1, 0, 0, -2, 0), dtype=object).reshape(7,7))
    Cycle2 = Quiver(2*np.array((0, 1, 0, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, -1, 0), dtype=object).reshape((7,7)))
    runBenchmarkMatchSets(list(generateMulticlassDataset([A7, X7, Cycle2], 2)), [A7, X7, Cycle2], useSurface=False)
    

if __name__ == "__main__":
    main()