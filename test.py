from Quiver import *
import unittest
import random

class QuiverInitializationTestCase(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver([[0,2], [-2,0]])
        self.quiver = self.quiver.mutate(1)

    def testLen(self):
        assert self.quiver.n == 2, 'incorrect rank n'
    
    def testVerts(self):
        assert self.quiver.vertices == [0,1], 'incorrect list of vertices'

    def testEdges(self):
        assert self.quiver.numEdges == 2, 'incorrect number of edges'

    def testDet(self):
        assert self.quiver.determinant() == 4, 'incorrect determinant'

    def testConnected(self):
        assert self.quiver.connected(), 'incorrect connectedness'

    def testAbundant(self):
        assert self.quiver.abundant(), 'incorrect abundance'


class QuiverMutation3vertTestCase(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver([[0,2,0], [-2,0,1], [0,-1,0]])
    
    def testInvolution(self):
        for i in self.quiver.vertices:
            assert self.quiver == self.quiver.mutate(i).mutate(i), 'mutation is not an involution'

    def testAcyclicCycle(self):
        assert self.quiver == self.quiver.mutate(0).mutate(1).mutate(2), 'source mutation cycle did not cycle'
        assert self.quiver == self.quiver.mutate(2).mutate(1).mutate(0), 'sink mutation cycle did not cycle'
    
    def testCommutingMutations(self):
        assert self.quiver == self.quiver.mutate(2).mutate(0).mutate(2).mutate(0), 'commuting mutation cycle did not cycle'
    
    def testA2Cycle(self):
        assert self.quiver.mutate(1).mutate(2).mutate(1).mutate(2).mutate(1) == self.quiver.mutate(2).mutate(1).mutate(2).mutate(1).mutate(2), 'A2 mutation cycle did not cycle'

    def testMutationInverseRandom(self):
        r = []
        for i in range(12):
            r.append(random.randrange(3))
        q = self.quiver
        for m in r:
            q = q.mutate(m)

        assert q.connected(), f'mutating at {r} resulted in a disconnected quiver'
        assert q.n == 3, f'mutating at {r} changed the number of vertices'

        for m in r[::-1]:
            q = q.mutate(m)
        assert self.quiver == q, f'applying the mutation sequence {r + r[::-1]} did not fix the quiver'

    def testCyclic(self):
        assert self.quiver.acyclic(), 'failed to identify acyclic'
        assert len(self.quiver.sources()) == 1 and len(self.quiver.sinks())==1, 'failed to identify sinks or sources'
        assert not self.quiver.mutate(1).acyclic(), 'failed to identify cyclic'
        assert len(self.quiver.mutate(1).sources()) == 0 and len(self.quiver.mutate(1).sinks())==0, f'incorrectly found sinks {self.quiver.mutate(1).sinks()} or sources {self.quiver.mutate(1).sources()}'

    def testOpposite(self):
        q = self.quiver.mutate(0).mutate(2)
        assert q == self.quiver.oppositeQuiver(), 'opposite quiver not generated'

    def testFork(self):
        q = self.quiver.mutate(1).mutate(0)
        assert q.forkWithPOR(0), 'failed to recognize fork'
        assert q.mutate(2).forkWithPOR(2), 'failed to recognize fork'



class QuiverMutation4vertTestCase(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver([[0,2,0,6], [-2,0,1,0], [0,-1,0,7], [-6,0,-7,0]])
    
    def testCommutingMutations(self):
        assert self.quiver == self.quiver.mutate(2).mutate(0).mutate(2).mutate(0), 'commuting mutation cycle did not cycle'
    
    def testA2Cycle(self):
        assert self.quiver.mutate(1).mutate(2).mutate(1).mutate(2).mutate(1) == self.quiver.mutate(2).mutate(1).mutate(2).mutate(1).mutate(2), 'A2 mutation cycle did not cycle'

    def testMutationInverseRandom(self):
        r = []
        for i in range(12):
            r.append(random.randrange(3))
        q = self.quiver
        for m in r:
            q = q.mutate(m)
        
        assert q.cyclic_order() != False, f'{r} mutated into something without a good cyclic ordering'
        assert q.determinant() == self.quiver.determinant(), f'determinant changed after mutating at {r}'
        assert not q.vortex(), f'acyclic quiver became a vortex after mutating at {r}'

        for m in r[::-1]:
            q = q.mutate(m)
        assert self.quiver == q, f'applying the mutation sequence {r + r[::-1]} did not fix the quiver'

    def testDet(self):
        for i in self.quiver.vertices:
            assert self.quiver.determinant() == self.quiver.mutate(i).determinant(), f'applying a mutation changed the determinant'

    def testPrefork(self):
        M = [1,0]
        q = self.quiver
        for m in M:
            q = q.mutate(m)
        assert q.forkWithPOR(0), f'applying mutation sequence {M} did not produce a prefork with por {M[:-1]}'
        
    def testFork(self):
        M = [1,0,3,2,0,1,2,3]
        q = self.quiver
        for m in M:
            q = q.mutate(m)
        assert q.forkWithPOR(3), f'applying mutation sequence {M} did not produce a fork with por {M[:-1]}'


class QuiverInvariants3dcTestCase(unittest.TestCase):
    def setUp(self):
        #These are disconnected quivers with only oriented cycles on 0,2,1.
        self.quiver_c1 = Quiver([[0,-3,6,0], [3,0,-5,0], [-6,5,0,0], [0,0,0,0]])
        self.quiver_c2 = Quiver([[0,3,-6,0], [-3,0,5,0], [6,-5,0,0], [0,0,0,0]])
        self.quiver_a1 = Quiver([[0,-2,6,0], [2,0,-5,0], [-6,5,0,0], [0,0,0,0]])
        self.quiver_a2 = Quiver([[0,2,-6,0], [-2,0,5,0], [6,-5,0,0], [0,0,0,0]])
        self.quivers = [self.quiver_c1, self.quiver_c2, self.quiver_a1, self.quiver_a2]
        self.areMutCyclic = [True, True, False, False]
        self.orientedCycles = [(0,2,1), (0,1,2), (0,2,1), (0,1,2)]
        self.markovs = [3**2+5**2+6**2 - 3*5*6, 3**2+5**2+6**2 - 3*5*6, 2**2+5**2+6**2 - 2*5*6, 2**2+5**2+6**2 - 2*5*6]

    def testAcyclic(self):
        for q in self.quivers:
            assert not q.acyclic()

    def testConnected(self):
        for q in self.quivers:
            assert not q.connected()

    def testCyclicSubquiver(self):
        for i in range(len(self.quivers)):
            if self.areMutCyclic[i]:
                assert self.quivers[i].hasMutCyclicSubquiver(), 'mutation cyclic subquiver is incorrectly unnoticed'
            else:
                assert not self.quivers[i].hasMutCyclicSubquiver(), 'mutation acyclic quiver is incorrectly marked as mutation cyclic'
        
    def testCycles(self):
        for q in self.quivers:
            assert [0,1,2] in list(q.chordless_cycles()) or [0,2,1] in list(q.chordless_cycles()), 'missing chordless cycle'
            assert len(list(q.chordless_cycles()))==1, 'incorrectly found more cycles'

    def testCyclicOrdering(self):
        sigmas = [q.cyclic_order() for q in self.quivers]
        assert False not in sigmas, "incorrectly failed to find a cyclic ordering"
        for i in range(len(self.quivers)):
            sigma = sigmas[i]
            c = self.orientedCycles[i]
            order_respected = [sigma[c[i]] < sigma[c[(i+1)%len(c)]] for i in range(len(c))]
            assert order_respected.count(False)==1, f"incorrect winding number ({order_respected.count(False)}) from cyclic ordering ({sigma}) for quiver {i}"

    def testMarkov(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].markov() == self.markovs[i], f"incorrect markov invariant {i, self.quivers[i].markov(), self.markovs[i]}"


class QuiverInvariants4AlexTestCase(unittest.TestCase):
    def setUp(self):
        self.cycleQuiver = lambda w : Quiver([[0,w[0],0,-w[3]],[-w[0],0,w[1],0],[0,-w[1],0,w[2]],[w[3],0,-w[2],0]])
        self.pathQuiver = lambda w : Quiver([[0,w[0],0,0],[-w[0],0,w[1],0],[0,-w[1],0,w[2]],[0,0,-w[2],0]])
        self.clawQuiver = lambda w : Quiver([[0,w[0],w[1],w[2]],[-w[0],0,0,0],[-w[1],0,0,0],[-w[2],0,0,0]])

    def testHelperFuns(self):
        assert self.pathQuiver([1,2,3]).mutate(0).mutate(1).mutate(2).mutate(3) == self.pathQuiver([1,2,3])
        assert self.pathQuiver([14,22,13]) == self.cycleQuiver([14,22,13,0])
        assert self.pathQuiver([4,0,0]) == self.clawQuiver([4,0,0])
    
    def testMarkovCycle(self):
        assert self.cycleQuiver([1,1,1,1]).markov() == 4 - 1
        assert self.cycleQuiver([1,2,1,1]).markov() == 4 + 3 - 2
        m = lambda w : sum([x**2 for x in w]) - w[0]*w[1]*w[2]*w[3]
        for w in [[1,1,2,1], [3,4,2,1], [9,9,9,9], [0,1,1,1], [7,3,2,-1], [100,1000,100000,10000000]]:
            assert m(w) == self.cycleQuiver(w).markov(), f"incorrect cycle quiver markov with weights {w}"

    def testMarkovPath(self):
        assert self.pathQuiver([1,1,1]).markov() == 3
        m = lambda w : sum([x**2 for x in w])
        for w in [[1,1,2], [3,2,1], [9,9,9], [0,1,1], [7,3,-11], [1000,100000,10000000]]:
            assert m(w) == self.pathQuiver(w).markov(), f"incorrect path quiver markov with weights {w}"

    def testMarkovClaw(self):
        assert self.clawQuiver([1,1,1]).markov() == 3
        m = lambda w : sum([x**2 for x in w])
        for w in [[1,1,2], [3,2,1], [9,9,9], [10,0,1], [7,3,-11], [1000,100000,10000000]]:
            assert m(w) == self.clawQuiver(w).markov(), f"incorrect claw quiver markov with weights {w}"

    def testAlexPath(self):
        assert self.pathQuiver([2,2,2]).alexander_poly() == polynomial.polynomial([1,8,-2,8,1])
        a = lambda w : polynomial.polynomial([1,-4,6,-4,1]) + (w[0]**2 + w[1]**2 + w[2]**2) * polynomial.polynomial([0,1,-2,1]) + (w[0]*w[2])**2 * polynomial.polynomial([0,0,1,0,0])
        for w in [[3,2,1], [100,1,99], [-19,-32,-100], [25,0,25]]:
            assert a(w) == self.pathQuiver(w).alexander_poly(), f"incorrect path quiver alexander poly with weights {w}"

    def testAlexClaw(self):
        assert self.clawQuiver([2,2,2]).alexander_poly() == polynomial.polynomial([1,8,-18,8,1])
        a = lambda w : polynomial.polynomial([1,-4,6,-4,1]) + (w[0]**2 + w[1]**2 + w[2]**2) * polynomial.polynomial([0,1,-2,1])
        for w in [[3,2,1], [100,1,99], [-19,-32,-100], [25,0,25]]:
            assert a(w) == self.clawQuiver(w).alexander_poly(), f"incorrect path quiver alexander poly with weights {w}"

    def testAlexCycle(self):
        assert self.cycleQuiver([3,3,3,3]).alexander_poly() == polynomial.polynomial([1,-4+36-81, 6+2*45, -4+36-81,1])
        a = lambda w: polynomial.polynomial([1,-4,6,-4,1]) + polynomial.polynomial([0,1,-2,1])*(sum([x**2 for x in w]) - w[0]*w[1]*w[2]*w[3]) + polynomial.polynomial([0,0,1])*(w[0]**2 *w[2]**2 + w[1]**2 * w[3]**2 - 2*w[0]*w[1]*w[2]*w[3])
        for w in [[1,1,1,1], [4,3,2,1], [9,1,1,9], [10000,3,10000,2], [2,2,2,4],[0,0,0,0], [19,23,41,-11]]:
            assert a(w) == self.cycleQuiver(w).alexander_poly(), f"incorrect cycle quiver alexander poly with weights {w}"


class QuiverInvariants3TestCase(unittest.TestCase):
    def setUp(self):
        #These are disconnected quivers with only oriented cycles on 0,2,1.
        self.quiver_c1 = Quiver([[0,-3,6], [3,0,-5], [-6,5,0]])
        self.quiver_c2 = Quiver([[0,-5,6], [5,0,-5], [-6,5,0]])
        self.quiver_a1 = Quiver([[0,2,-6], [-2,0,5], [6,-5,0]])
        self.quiver_a2 = Quiver([[0,-3,6], [3,0,-18], [-6,18,0]])
        self.quivers = [self.quiver_c1, self.quiver_c2, self.quiver_a1, self.quiver_a2]
        self.areMutCyclic = [True, True, False, False]
        self.orientedCycles = [(0,2,1), (0,2,1), (0,1,2), (0,2,1)]
        self.markovs = [3**2+5**2+6**2 - 3*5*6, 5**2+5**2+6**2 - 5*5*6, 2**2+5**2+6**2 - 2*5*6, 3**2 + 6**2]

    def testAcyclic(self):
        for q in self.quivers:
            assert not q.acyclic()

    def testConnected(self):
        for q in self.quivers:
            assert q.connected()

    def testCycles(self):
        for q in self.quivers:
            assert [0,1,2] in list(q.chordless_cycles()) or [0,2,1] in list(q.chordless_cycles()), 'missing chordless cycle'
            assert len(list(q.chordless_cycles()))==1, 'incorrectly found more cycles'

    def testCyclicOrdering(self):
        sigmas = [q.cyclic_order() for q in self.quivers]
        assert False not in sigmas, "incorrectly failed to find a cyclic ordering"
        for i in range(len(self.quivers)):
            sigma = sigmas[i]
            c = self.orientedCycles[i]
            order_respected = [sigma[c[i]] < sigma[c[(i+1)%len(c)]] for i in range(len(c))]
            assert order_respected.count(False)==1, f"incorrect winding number ({order_respected.count(False)}) from cyclic ordering ({sigma}) for quiver {i}"

    def testMarkov(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].markov() == self.markovs[i], f"incorrect markov invariant {i, self.quivers[i].markov(), self.markovs[i]}"

    def testAlexander(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].alexander_poly() == polynomial.polynomial([-1, 3-self.markovs[i], -(3-self.markovs[i]),1]), f"incorrect alexander poly {str(self.quivers[i].alexander_poly()), self.markovs[i]}"

class RealizableAlexandersTestCase(unittest.TestCase):
    def testSearchA1(self):
        list_of_A1 = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,1])))
        assert len(list_of_A1)==1, "acyclics from alexander poly is wrong size"
        assert list_of_A1[0].matrix == np.array([0], dtype='object')

    def testSearchDisjoint4(self):
        list_of_A1_A1_A1_A1 = list(generate_acyclics_from_Alexander(polynomial.polynomial([1,-4,6,-4,1,0,0])))
        assert len(list_of_A1_A1_A1_A1)==1, "acyclics from alexander poly is wrong size"
        assert not list_of_A1_A1_A1_A1[0].connected(), "incorrectly found connected quiver"
        assert list_of_A1_A1_A1_A1[0].matrix == 0 * list_of_A1_A1_A1_A1[0].matrix, "incorrect acyclic found"

    def testSearchSmall3(self):
        disconn = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,2,-2,1])))
        assert all([not x.connected for x in disconn])
        assert all([x.alexander_poly() == polynomial.polynomial([-1,2,-2,1]) for x in disconn])

    def testSearchMid3(self):
        conn = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,-17,17,1])))
        assert all([x.connected for x in conn])
        assert all([x.alexander_poly() == polynomial.polynomial([-1,-17,17,1]) for x in conn])
        assert len(conn) == 4 #wts: 321, 312, 420, 222

    

class QuiverHashing(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver([[0,-3,6,8], [3,0,-5,9], [-6,5,0,-11], [-8,-9,11,0]])
        self.quiver_dup = copy.deepcopy(self.quiver)

    def testHash(self):
        assert hash(self.quiver) == hash(self.quiver_dup), "hashing cares about data besides the weights"
        assert hash(self.quiver.mutate(1)) == hash(self.quiver_dup.mutate(1)), "hashing cares about data besides the weights"

    def testMutChanges(self):
        for i in self.quiver.vertices:
            assert hash(self.quiver.mutate(i)) != hash(self.quiver), f"mutation at {i} fixed the hash"
        

class PermutationTests(unittest.TestCase):
    def setUp(self):
        self.s3 = [tuple(p) for p in permutations(3)]
        self.s2 = [tuple(p) for p in permutations(2)]
    
    def testSizes(self):
        assert len(self.s3)==6
        assert len(self.s2)==2
        assert len(set(self.s3))==6
        assert len(set(self.s2))==2
        assert len(self.s3[0])==3
        assert len(self.s2[0])==2

    def testMembers(self):
        assert (0,1) in self.s2
        assert (1,0) in self.s2
        assert (1,2,0) in self.s3
        assert (0,2,1) in self.s3

    def testParity(self):
        assert sign_perm([0,1]) == 1
        assert sign_perm([1,0]) == -1
        assert sign_perm([1,2,0]) == 1
        assert sign_perm([0,2,1]) == -1


if __name__ == "__main__":
    unittest.main()