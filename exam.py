from quiver import *
from prediction import *
from generation import *

def main():
    """ 
    Evaluation of several previous quiver mutation databases using the baselines in prediction.
     """
    # ML for AG: A7, X7, and the oriented cycle with weights = 2
    #  standard exchange matrices were imported using sagemath vector() on the associated B-matrix.
    A7 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0), dtype=object).reshape((7,7)))
    X7 = Quiver(np.array((0, 2, -1, 0, 0, 0, 0, -2, 0, 1, 0, 0, 0, 0, 1, -1, 0, 1, -1, 1, -1, 0, 0, -1, 0, 2, 0, 0, 0, 0, 1, -2, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 1, 0, 0, -2, 0), dtype=object).reshape(7,7))
    Cycle2 = Quiver(2*np.array((0, 1, 0, 0, 0, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 1, 0, 0, 0, 0, -1, 0), dtype=object).reshape((7,7)))
    print("baseline for A7, X7, and Cycle2")
    runBaselineMatchSets(list(generateMulticlassDataset([A7, X7, Cycle2], 2)), [A7, X7, Cycle2], useSurface=False)
    # Results: 
    # Ignoring Alexander Polynomial for non-totally-proper reps: 602 (66.66666666666667%)
    # ... and uniform guessing: 752.5 (83.33333333333333%)
    # Using Alexander Polynomials for the cyclic examples: 903 (100.0%)
    # ... and uniform guessing: 903.0 (100.0%)
    # they reported: 55%

    # QM, SD, and ML: A7, X7, and Affine E6
    AffineE6 = Quiver(np.array((0, 1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0), dtype=object).reshape(7,7))
    print("baseline for A7, X7, and Affine E6")
    runBaselineMatchSets(list(generateMulticlassDataset([A7, X7, AffineE6], 2)), [A7, X7, AffineE6], useSurface=False)
    # Results:
    # Ignoring Alexander Polynomial for non-totally-proper reps: 301 (33.333333333333336%)
    # ... and uniform guessing: 602.0 (66.66666666666667%)
    # Using Alexander Polynomials for the cyclic examples: 903 (100.0%)
    # ... and uniform guessing: 903.0 (100.0%)
    # they reported: 31%

    # MLMA: A4, A4Doubled, A4doubleMiddle, D4, D4DoubleOne, Box, Vortex, (Figures A1, A2, A4, A5, A6)
    A4 = A7.subquiver([0,1,2,3])
    A4Doubled = Quiver(2*A4.matrix)
    A4DoubleMiddle = Quiver([[0,1,0,0],[-1,0,-2,0],[0,2,0,1],[0,0,-1,0]])
    D4 = Quiver([[0,1,1,1],[-1,0,0,0],[-1,0,0,0],[-1,0,0,0]])
    D4DoubleOne = Quiver([[0,1,1,2],[-1,0,0,0],[-1,0,0,0],[-2,0,0,0]])
    box = Quiver([[0,2,0,-2],[-2,0,2,0],[0,-2,0,2],[2,0,-2,0]])
    lollipop = Quiver([[0,2,-2,0],[-2,0,2,2],[2,-2,0,0],[0,-2,0,0]])
    print("baseline for variations on A4, D4, as well as the box quiver and lollipop.")
    runBaselineMatchSets(list(generateMulticlassDataset([A4, A4Doubled, A4DoubleMiddle, D4, D4DoubleOne, box, lollipop], 5)), [A4, A4Doubled, A4DoubleMiddle, D4, D4DoubleOne, box, lollipop], useSurface=False)
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
    print("baseline for A,D,E,DE 11")
    runBaselineMatchSets(list(generateMulticlassDataset([A11, D11, E11, DE11], 1)), [A11, D11, E11, DE11], useSurface=False)
    # Ignoring Alexander Polynomial for non-totally-proper reps: 484 (100.0%)
    # ... and uniform guessing: 484.0 (100.0%)  #(All acyclic, so all use Alexander)
    # Using Alexander Polynomials for the cyclic examples: 484 (100.0%)
    # ... and uniform guessing: 484.0 (100.0%)
    # Reported: 93%

    # MLMA: classifying Mutation acyclic or not for 4-vertex quivers with small weights
    print("baseline all small weight 4-vertex quivers for Mutation-acyclicity")
    AllSmall4VertexData = [(r["quiver exchange matrix"], r["mutation acyclic"]) for r in readCSVDatabase("dataset/Labeled_4vertex_2MaxWeight_depth5.csv")]
    runBaselineAcyclicLocal(AllSmall4VertexData)
    #Local tests
    # data balance (a vs c): 123562 vs 106846
    # Unknown->incorrect: 57758 (25.067705982431164%)
    # *Unknown->True: 172328 (74.79254192562759%)
    # Unknown->False: 115838 (50.27516405680358%)

    # The same, using some training data from the dataset at random
    runBaselineAcyclicWithTraining(AllSmall4VertexData)
    # Comparing against a train set.
    # Found 124 distinct cyclic invariants and 204 distinct acyclic invariants.
    # They intersect for 1 values.

    # Ignoring Alexander Polynomial for cyclic examples:
    # Unknown->incorrect: 230194 (99.90712128051109%)
    # *Unknown->True: 230378 (99.98697961876324%)
    # Unknown->False: 230224 (99.92014166174786%)
    # Using Alexander Polynomials for the cyclic examples:
    # Unknown->incorrect: 230194 (99.90712128051109%)
    # Unknown->True: 230378 (99.98697961876324%)
    # Unknown->False: 230224 (99.92014166174786%)

    def TypeAQuiver(n):
        B = np.zeros((n,n), dtype=object)
        for i in range(n-1):
            B[i,i+1] = 1
            B[i+1,i] = -1
        return Quiver(B)

    def TypeDQuiver(n):
        assert n>=4
        B = np.zeros((n,n), dtype=object)
        B[:3,:3] = [[0,0,1],[0,0,1],[-1,-1,0]]
        B[2:,2:] = TypeAQuiver(n-2).matrix
        return Quiver(B)

    def TypeEQuiver(n):
        assert n>=6
        B = np.zeros((n,n), dtype=object)
        B[:5,:5] = [[0,1,0,0,0],[-1,0,0,0,1],[0,0,0,1,0],[0,0,-1,0,1],[0,-1,0,-1,0]]
        B[4:,4:] = TypeAQuiver(n-4).matrix
        return Quiver(B)
        
    # He et. al. 2025
    print("baseline for A,D,E for 3<=n<9")
    ManyTrees = [TypeAQuiver(n) for n in range(3,9)] + [TypeDQuiver(n) for n in range(4,9)] + [TypeEQuiver(n) for n in range(6,9)]
    runBaselineMatchSets(list(generateMulticlassDataset(ManyTrees, 2)), ManyTrees, useSurface=False)
    # not run for n>8 or more mutations, as the alexander polynomial distinguishes all of these. See the thesis of Amanda Shwartz for a proof.
    # Ignoring Alexander Polynomial for non-totally-proper reps: 3164 (100.0%)
    # ... and uniform guessing: 3164.0 (100.0%)
    # Using Alexander Polynomials for the cyclic examples: 3164 (100.0%)
    # ... and uniform guessing: 3164.0 (100.0%)
    


if __name__ == "__main__":
    main()