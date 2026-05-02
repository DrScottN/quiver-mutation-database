from quiver import *
from prediction import *
from generation import *

def main():
    """ 
    Evaluation of several previous quiver mutation databases using the benchmarks in prediction.
     """
    # ML for AG: A7, X7, and the oriented cycle with weights = 2
    #  standard exchange matrices were imported using sagemath vector() on the associated B-matrix.
    A7 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0), dtype=object).reshape((7,7)))
    X7 = Quiver(np.array((0, 2, -1, 0, 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, -1, 0, 1, -1, 1, -1, 0, 0, -1, 0, 2, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 1, 0, 0, -2, 0), dtype=object).reshape(7,7))
    Cycle2 = Quiver(2*np.array((0, 1, 0, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, -1, 0), dtype=object).reshape((7,7)))
    #runBenchmarkMatchSets(list(generateMulticlassDataset([A7, X7, Cycle2], 2)), [A7, X7, Cycle2], useSurface=False)
    # Results: 
    # Ignoring Alexander Polynomial for non-totally-proper reps: 602 (66.66666666666667%)
    # ... and uniform guessing: 752.5 (83.33333333333333%)
    # Using Alexander Polynomials for the cyclic examples: 903 (100.0%)
    # ... and uniform guessing: 903.0 (100.0%)
    # they reported: 55%

    # QM, SD, and ML: A7, X7, and Affine E6
    AffineE6 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0), dtype=object).reshape(7,7))
    # runBenchmarkMatchSets(list(generateMulticlassDataset([A7, X7, AffineE6], 2)), [A7, X7, AffineE6], useSurface=False)
    # Results:
    # Ignoring Alexander Polynomial for non-totally-proper reps: 301 (33.333333333333336%)
    # ... and uniform guessing: 602.0 (66.66666666666667%)
    # Using Alexander Polynomials for the cyclic examples: 903 (100.0%)
    # ... and uniform guessing: 903.0 (100.0%)
    # they reported: 31%

    # QM, SD, and ML: (3,3,3) cycle, (2,2,4) cycle (acyclic), and a collection of random 3x3 skew-symm matrices.
    #  skipped, they don't say what 'random' means. Note that markov+gcd distinguishes almost all of these.

    # MLMA: A4, A4Doubled, A4doubleMiddle, D4, D4DoubleOne, Box, Vortex, (Figures A1, A2, A4, A5, A6)
    A4 = A7.subquiver([0,1,2,3])
    A4Doubled = Quiver(2*A4.matrix)
    A4DoubleMiddle = Quiver([[0,1,0,0],[-1,0,-2,0],[0,2,0,1],[0,0,-1,0]])
    D4 = Quiver([[0,1,1,1],[-1,0,0,0],[-1,0,0,0],[-1,0,0,0]])
    D4DoubleOne = Quiver([[0,1,1,2],[-1,0,0,0],[-1,0,0,0],[-2,0,0,0]])
    box = Quiver([[0,2,0,-2],[-2,0,2,0],[0,-2,0,2],[2,0,-2,0]])
    lollipop = Quiver([[0,2,-2,0],[-2,0,2,2],[2,-2,0,0],[0,-2,0,0]])
    #runBenchmarkMatchSets(list(generateMulticlassDataset([A4, A4Doubled, A4DoubleMiddle, D4, D4DoubleOne, box, lollipop], 5)), [A4, A4Doubled, A4DoubleMiddle, D4, D4DoubleOne, box, lollipop], useSurface=False)
    # Results:
    # Ignoring Alexander Polynomial for non-totally-proper reps: 8736 (85.71428571428571%)
    # ... and uniform guessing: 9464.0 (92.85714285714286%)
    # Using Alexander Polynomials for the cyclic examples: 10192 (100.0%)
    # ... and uniform guessing: 10192.0 (100.0%)
    # Reported: 61%

    #ML meets AC: A,D,E,DE on 11 vertices
    # *Note that they consider skew-symmetrizable, and so have a larger dataset. We expect our results easily generalize.
    A11 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)).reshape((11,11)))
    D11 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0)).reshape((11,11)))
    E11 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)).reshape((11,11)))
    DE11 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)).reshape((11,11)))
    # runBenchmarkMatchSets(list(generateMulticlassDataset([A11, D11, E11, DE11], 1)), [A11, D11, E11, DE11], useSurface=False)
    # Ignoring Alexander Polynomial for non-totally-proper reps: 484 (100.0%)
    # ... and uniform guessing: 484.0 (100.0%)  #(All acyclic, so all use Alexander)
    # Using Alexander Polynomials for the cyclic examples: 484 (100.0%)
    # ... and uniform guessing: 484.0 (100.0%)
    # Reported: 93%

    # MLMA: classifying Mutation acyclic or not for 4-vertex quivers with small weights
    #todo

if __name__ == "__main__":
    main()