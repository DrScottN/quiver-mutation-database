from Quiver import *
from unknown import Unknown
import unittest
import random
import itertools

class QuiverInitializationTestCase(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver([[0,2], [-2,0]])
        self.quiver = self.quiver.mutate(1)
        self.quiver_iso = Quiver([[0,0],[0,0]])

    def testLen(self):
        assert self.quiver.n == 2, 'incorrect rank n'
        assert self.quiver_iso.n == 2, 'incorrect rank n'
    
    def testVerts(self):
        assert self.quiver.vertices == [0,1], 'incorrect list of vertices'
        assert self.quiver_iso.vertices == [0,1], 'incorrect list of vertices'

    def testEdges(self):
        assert self.quiver.numEdges == 2, 'incorrect number of edges'
        assert self.quiver_iso.numEdges == 0, 'incorrect number of edges'

    def testDet(self):
        assert self.quiver.determinant() == 4, 'incorrect determinant'
        assert self.quiver_iso.determinant() == 0, 'incorrect determinant'

    def testConnected(self):
        assert self.quiver.connected(), 'incorrect connectedness'
        assert not self.quiver_iso.connected(), 'incorrect connectedness'

    def testAbundant(self):
        assert self.quiver.abundant(), 'incorrect abundance'
        assert not self.quiver_iso.abundant(), 'incorrect abundance'

class QuiverConnectedIsolatedTestCase(unittest.TestCase):
    def testA1Isolated(self):
        q = Quiver([[0,1,0],[-1,0,0],[0,0,0]])
        assert not q.connected()

    def testmA1Isolated(self):
        q = Quiver([[0,-1,0],[1,0,0],[0,0,0]])
        assert not q.connected()



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
        self.incomplete_quiver = Quiver(np.matrix([[0,2,-3,0],[-2,0,2,0],[3,-2,0,0], [0,0,0,0]]))
        self.vortex_quiver = Quiver(np.matrix([[0,2,-4,-6],[-2,0,2,-6],[4,-2,0,-6], [6,6,6,0]]))
        self.A4 = Quiver([[0,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,0]])
    
    def testCommutingMutations(self):
        assert self.quiver == self.quiver.mutate(2).mutate(0).mutate(2).mutate(0), 'commuting mutation cycle did not cycle'
        assert self.incomplete_quiver.mutate(1) == self.incomplete_quiver.mutate(3).mutate(1), 'commuting mutation cycle did not cycle'
        assert self.incomplete_quiver.mutate(1) == self.incomplete_quiver.mutate(1).mutate(3), 'commuting mutation cycle did not cycle'
        assert self.vortex_quiver.mutate(1).mutate(2) != self.vortex_quiver.mutate(2).mutate(1), 'too many mutations commute'
        assert self.A4.mutate(0).mutate(3) == self.A4.mutate(3).mutate(0)

    def testA2Cycle(self):
        assert self.quiver.mutate(1).mutate(2).mutate(1).mutate(2).mutate(1) == self.quiver.mutate(2).mutate(1).mutate(2).mutate(1).mutate(2), 'A2 mutation cycle did not cycle'
        assert self.A4.mutate(1).mutate(2).mutate(1).mutate(2).mutate(1) == self.A4.mutate(2).mutate(1).mutate(2).mutate(1).mutate(2)

    def testMutationInverseRandom(self):
        r = []
        for i in range(12):
            r.append(random.randrange(4))
        q = self.quiver
        for m in r:
            q = q.mutate(m)
        
        assert q.cyclic_order() != False, f'{r} mutated into something without a good cyclic ordering'
        assert q.determinant() == self.quiver.determinant(), f'determinant changed after mutating at {r}'
        assert not q.vortex(), f'acyclic quiver became a vortex after mutating at {r}'
        assert q.vortex_free(), f'vortex free disagrees with .vortex()'

        for m in r[::-1]:
            q = q.mutate(m)
        assert self.quiver == q, f'applying the mutation sequence {r + r[::-1]} did not fix the quiver'

    def testMutationInverseRandom_incomplete(self):
        r = []
        for i in range(12):
            r.append(random.randrange(4))
        q = self.incomplete_quiver
        for m in r:
            q = q.mutate(m)
        
        assert q.cyclic_order() != False, f'{r} mutated into something without a good cyclic ordering'
        assert q.determinant() == self.incomplete_quiver.determinant(), f'determinant changed after mutating at {r}'
        assert not q.vortex(), f'acyclic quiver became a vortex after mutating at {r}'
        assert q.vortex_free(), f'vortex free disagrees with .vortex()'

        for m in r[::-1]:
            q = q.mutate(m)
        assert self.incomplete_quiver == q, f'applying the mutation sequence {r + r[::-1]} did not fix the quiver'

    def testDet(self):
        for i in self.quiver.vertices:
            assert self.quiver.determinant() == self.quiver.mutate(i).determinant(), f'applying a mutation changed the determinant'
            assert self.incomplete_quiver.determinant() == self.incomplete_quiver.mutate(i).determinant(), f'applying a mutation changed the determinant'
            assert self.vortex_quiver.determinant() == self.vortex_quiver.mutate(i).determinant(), f'applying a mutation changed the determinant'
            assert self.A4.determinant() == self.A4.mutate(i).determinant()

    def testPrefork(self):
        M = [1,0]
        q = self.quiver
        for m in M:
            q = q.mutate(m)
        assert q.preForkWithPOR(M[-1]), f'applying mutation sequence {M} did not produce a prefork with por {M[-1]}'
        q = self.incomplete_quiver
        for m in M:
            q = q.mutate(m)
        assert not any([q.preForkWithPOR(i) for i in range(q.n)]), f'incorrectly marked a disconnected quiver as a prefork'

    def testPrefork2(self):
        M = [1,0,1]
        q = self.quiver
        for m in M:
            q = q.mutate(m)
        assert q.preForkWithPOR(M[-1]), f'applying mutation sequence {M} did not produce a prefork with por {M[-1]}'
        q = self.incomplete_quiver
        for m in M:
            q = q.mutate(m)
        assert not any([q.preForkWithPOR(i) for i in range(q.n)]), f'incorrectly marked a disconnected quiver as a prefork'

    def testPrefork3(self):
        q = self.incomplete_quiver
        q = q.mutate(2)
        for i,j in itertools.combinations(q.vertices, 2):
            assert not q.preForkWithVertices(2, i, j), f'incorrected found a prefork with vertices {i,j}'
            assert q.preForkWithVertices(2, i, j) == q.preForkWithVertices(2, j, i), f'incorrectly cares about order of prefork vertices {i,j}'
        assert not q.preForkWithPOR(2)
        
    def testFork(self):
        M = [1,0,3,2,0,1,2,3]
        q = self.quiver
        for m in M:
            q = q.mutate(m)
        assert q.forkWithPOR(M[-1]), f'applying mutation sequence {M} did not produce a fork with por {M[-1]}'
        q = self.incomplete_quiver
        for m in M:
            q = q.mutate(m)
        assert not any([q.forkWithPOR(i) for i in range(q.n)]), f'incorrectly marked a disconnected quiver as a fork'
        q = self.vortex_quiver
        for m in M:
            q = q.mutate(m)
        assert q.forkWithPOR(M[-1]), f'applying mutation sequence {M} did not produce a fork with por {M[-1]}'

    def testFork2(self):
        assert not self.incomplete_quiver.forkWithPOR(1)
        assert not self.quiver.forkWithPOR(2)
        assert not any([self.vortex_quiver.forkWithPOR(i) for i in self.vortex_quiver.vertices])
        assert not any([self.vortex_quiver.mutate(3).forkWithPOR(i) for i in self.vortex_quiver.vertices])
        assert not any([self.vortex_quiver.mutate(1).forkWithPOR(i) for i in self.vortex_quiver.vertices])
        assert self.vortex_quiver.mutate(0).forkWithPOR(0) 
        assert self.vortex_quiver.mutate(2).forkWithPOR(2) 

    def testChordlessCycles(self):
        assert set(self.A4.chordless_cycles()) == set()
        assert set(self.vortex_quiver.chordless_cycles()) == set([(0,1,2),(0,1,3),(0,2,3),(1,2,3)]), f"missing chordless cycles {set(self.vortex_quiver.chordless_cycles())}"

    def testWindingDictionary(self):
        assert self.vortex_quiver.winding_data([3,2,1,0])[(0,1,2)] == (2,0), f"incorrect winding, got {self.vortex_quiver.winding_data([3,2,1,0])[(0,1,2)]}"
        assert self.vortex_quiver.winding_data([3,2,1,0])[(0,1,3)] == (1,1), f"incorrect winding, got {self.vortex_quiver.winding_data([3,2,1,0])[(0,1,3)]}"
        assert self.vortex_quiver.winding_data([3,2,1,0])[(0,2,3)] == (0,2), f"incorrect winding, got {self.vortex_quiver.winding_data([3,2,1,0])[(0,2,3)]}"
        assert self.vortex_quiver.winding_data([3,2,1,0])[(1,2,3)] == (1,1), f"incorrect winding, got {self.vortex_quiver.winding_data([3,2,1,0])[(1,2,3)]}"

    def testProper(self):
        assert self.incomplete_quiver.proper()
        assert self.quiver.proper()
        assert not self.vortex_quiver.proper(), f"incorrectly found proper ordering {self.vortex_quiver.proper()} of a vortex"

    def testCasals(self):
        assert self.quiver.casals_det()==0
        assert self.incomplete_quiver.casals_det()==0
        assert self.vortex_quiver.casals_det()==0
        assert self.A4.casals_det()==1
        for i in range(self.quiver.n):
            assert self.vortex_quiver.mutate(i).casals_det()==0
            assert self.A4.mutate(i).casals_det()==1

    def testSeven(self):
        for i in range(self.quiver.n):
            assert self.A4.seven_congruence()==self.A4.mutate(i).seven_congruence()
            assert self.quiver.seven_congruence()==self.quiver.mutate(i).seven_congruence()
            assert self.vortex_quiver.seven_congruence()==self.vortex_quiver.mutate(i).seven_congruence()

    def testVortex(self):
        assert self.A4.vortex_free()
        assert not self.A4.vortex()
        assert not self.vortex_quiver.vortex_free()
        assert self.vortex_quiver.mutate(2).vortex_free()
        assert self.incomplete_quiver.vortex_free()
        assert not self.incomplete_quiver.vortex()

class ForkEdgeTestCase(unittest.TestCase):
    def setUp(self):
        self.iso2 = Quiver(np.matrix([[0,0],[0,0]]))
        self.twoForks = Quiver(np.matrix([[0,5,-15,0,0,0],[-5,0,2,0,0,0],[15,-2,0,0,0,0],[0,0,0,0,6,-105],[0,0,0,-6,0,10],[0,0,0,105,-10,0]]))
        self.prefork = Quiver(np.matrix([[0,1,4,-17],[-1,0,5,-18],[-4,-5,0,3],[17,18,-3,0]]))
        self.forkWithExtra = Quiver([[0,-3,6,0], [3,0,-5,0], [-6,5,0,0], [0,0,0,0]])

    def testNotForks(self):
        assert all([not self.iso2.forkWithPOR(i) for i in range(self.iso2.n)])
        assert all([not self.twoForks.forkWithPOR(i) for i in range(self.twoForks.n)])
        assert all([not self.prefork.forkWithPOR(i) for i in range(self.prefork.n)])
        assert all([not self.forkWithExtra.forkWithPOR(i) for i in range(self.forkWithExtra.n)])

    def testNotForksMutated(self):
        Q = self.iso2.mutate(1)
        assert all([not Q.forkWithPOR(i) for i in range(Q.n)])
        Q = self.twoForks.mutate(4).mutate(5).mutate(1)
        assert all([not Q.forkWithPOR(i) for i in range(Q.n)])
        Q = self.prefork.mutate(1).mutate(3).mutate(2)
        assert Q.forkWithPOR(2)
        Q = self.prefork.mutate(1).mutate(0)
        assert all([not Q.forkWithPOR(i) for i in range(Q.n)])

        

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
        self.gcd_vecs = [(3,1,1,0), (3,1,1,0), (2,1,1,0), (2,1,1,0)]#xx

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
            assert (0,1,2) in list(q.chordless_cycles()) or (0,2,1) in list(q.chordless_cycles()), 'missing chordless cycle'
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

    def testGCDvecs(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].gcd_vector() == self.gcd_vecs[i], f"incorrect gcd vector, {i, self.gcd_vecs[i], self.quivers[i].gcd_vector()}"


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


class SubquiverTestCase(unittest.TestCase):
    def setUp(self):
        self.markov_vortex = Quiver([[0,2,-2,4],[-2,0,2,7],[2,-2,0,9],[-4,-7,-9,0]])
        self.mutated_vortex = self.markov_vortex.mutate(1).mutate(3).mutate(2).mutate(0)
        self.small_weight = Quiver([[0,2,-1,-1],[-2,0,2,0],[1,-2,0,1],[1,0,-1,0]])
        self.mutated_small = self.small_weight.mutate(3)
        self.big_problem = Quiver([[0, -1900534033, 89268872382, -568349], [1900534033, 0, 17905678993253867017, -114000260996770], [-89268872382, -17905678993253867017, 0, 157067], [568349, 114000260996770, -157067, 0]] )
        self.markov_manual = lambda M : (M[0,1]**2) + (M[1,2]**2) + (M[2,0]**2) - abs(M[0,1]*M[1,2]*M[2,0])

    def testCyclicSubquivers(self):
        assert self.markov_vortex.hasMutCyclicSubquiver()
        assert not self.mutated_vortex.hasMutCyclicSubquiver()
        assert not self.small_weight.hasMutCyclicSubquiver()
        assert self.mutated_small.hasMutCyclicSubquiver()
        assert not self.big_problem.hasMutCyclicSubquiver()

    def testMarkov(self):
        for i in self.markov_vortex.vertices:
            R = self.markov_vortex.subquiverRemoveOneVertex(i)
            if not R.acyclic():
                assert self.markov_manual(R.matrix) == R.markov()
        
        for i in self.big_problem.vertices:
            R = self.big_problem.subquiverRemoveOneVertex(i)
            if not R.acyclic():
                assert self.markov_manual(R.matrix) == R.markov(), f"incorrect markov for {i, R.matrix}, {self.markov_manual(R.matrix)} vs {R.markov()}"
        
    def testAcyclicBig(self):
        assert self.big_problem.subquiverRemoveOneVertex(3).acyclic()
        assert self.big_problem.subquiverRemoveOneVertex(2).acyclic()
        assert not self.big_problem.subquiverRemoveOneVertex(1).acyclic()
        assert not self.big_problem.subquiverRemoveOneVertex(0).acyclic()

    def testAbundance(self):
        assert self.mutated_vortex.abundant()
        assert self.markov_vortex.abundant()
        assert not self.small_weight.abundant()
        assert not self.small_weight.abundant()
    

class QuiverInvariants3TestCase(unittest.TestCase):
    def setUp(self):
        self.quiver_c1 = Quiver([[0,-3,6], [3,0,-5], [-6,5,0]])
        self.quiver_c2 = Quiver([[0,-5,6], [5,0,-5], [-6,5,0]])
        self.quiver_a1 = Quiver([[0,2,-6], [-2,0,5], [6,-5,0]])
        self.quiver_a2 = Quiver([[0,-3,6], [3,0,-18], [-6,18,0]])
        self.quivers = [self.quiver_c1, self.quiver_c2, self.quiver_a1, self.quiver_a2]
        self.areMutCyclic = [True, True, False, False]
        self.orientedCycles = [(0,2,1), (0,2,1), (0,1,2), (0,2,1)]
        self.markovs = [3**2+5**2+6**2 - 3*5*6, 5**2+5**2+6**2 - 5*5*6, 2**2+5**2+6**2 - 2*5*6, 3**2 + 6**2]
        self.gcd_vecs = [(3,1,1), (1,5,1), (2,1,1), (3,3,6)]
        self.casals_dets = [0, 0, 2, 2]
        self.sevens_kers = [1, 1, 0, 0]

    def testAcyclic(self):
        for q in self.quivers:
            assert not q.acyclic()

    def testConnected(self):
        for q in self.quivers:
            assert q.connected()

    def testCycles(self):
        for q in self.quivers:
            assert (0,1,2) in list(q.chordless_cycles()) or (0,2,1) in list(q.chordless_cycles()), 'missing chordless cycle'
            assert len(list(q.chordless_cycles()))==1, 'incorrectly found more cycles'
        for _,d in self.quiver_a1.winding_data([2,1,0]).items():
            assert d in [(-2,3), (2,0)]


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

    def testGCDvecs(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].gcd_vector() == self.gcd_vecs[i], f"incorrect gcd vector, {i, self.gcd_vecs[i], self.quivers[i].gcd_vector()}"

    def testCasals(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].casals_det() == self.casals_dets[i]

    def testSeven(self):
        for i in range(len(self.quivers)):
            assert self.quivers[i].seven_congruence() == self.sevens_kers[i]

class RealizableAlexandersTestCase(unittest.TestCase):
    def testSearchA1(self):
        list_of_A1 = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,1])))
        assert len(list_of_A1)==1, "acyclics from alexander poly is wrong size"
        assert list_of_A1[0].matrix == np.array([0], dtype='object')

    def testSearchDisjoint4(self):
        list_of_A1_A1_A1_A1 = list(generate_acyclics_from_Alexander(polynomial.polynomial([1,-4,6,-4,1,0,0])))
        assert len(list_of_A1_A1_A1_A1)==1, "acyclics from alexander poly is wrong size"
        assert not list_of_A1_A1_A1_A1[0].connected(), "incorrectly found connected quiver"
        assert (list_of_A1_A1_A1_A1[0].matrix == 0 * list_of_A1_A1_A1_A1[0].matrix).all(), "incorrect acyclic found"

    def testSearchSmall3(self):
        disconn = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,2,-2,1])))
        assert all([x.alexander_poly() == polynomial.polynomial([-1,2,-2,1]) for x in disconn])
        assert all([not x.connected() for x in disconn])
        assert len(disconn) == 1, f"incorrect solution set {[Q.matrix for Q in disconn]}" #wts: 100


    def testSearchMid3(self):
        conn = list(generate_acyclics_from_Alexander(polynomial.polynomial([-1,-17,17,1])))
        assert all([x.connected() for x in conn])
        assert all([x.alexander_poly() == polynomial.polynomial([-1,-17,17,1]) for x in conn])
        assert len(conn) == 4, f"incorrect solution set {[Q.matrix for Q in conn]}" #wts: 321, 312, 420, 222

class QuiverIsoClass(unittest.TestCase):
    def setUp(self):
        self.quiver1 = Quiver([[0,3,2],[-3,0,1],[-2,-1,0]])
        self.quiver2 = Quiver([[0,3,1],[-3,0,2],[-1,-2,0]])
        self.Iso1 = isomorphismClass(self.quiver1, permutations(3))
        self.Iso2 = isomorphismClass(self.quiver2, permutations(3))

    def testDistinct(self):
        assert self.quiver1 in self.Iso1, "incorrectly misses quiver from its own iso class"
        assert self.quiver2 in self.Iso2, "incorrectly misses quiver from its own iso class"
        assert self.quiver2 not in self.Iso1, "incorrectly includes quiver from distinct iso class"
        assert self.quiver1 not in self.Iso2, "incorrectly includes quiver from distinct iso class"

    def testMarkov(self):
        assert self.quiver1.markov() == self.quiver2.markov()

    def testSinks(self):
        assert self.quiver1 in sink_set(self.quiver1)
        assert self.quiver2 in sink_set(self.quiver2)
        assert all([self.quiver2 not in isomorphismClass(R, permutations(3)) for R in sink_set(self.quiver1)])
        assert all([self.quiver1 not in isomorphismClass(R, permutations(3)) for R in sink_set(self.quiver2)])

    def testMinRep(self):
        B = -self.quiver1.matrix
        Qp = Quiver(B)
        for R in sink_set(self.quiver2):
            assert not any([Rp < Qp for Rp in isomorphismClass(R, permutations(3)) if (np.tril(Rp.matrix) >= 0).all()])

    def testMinRep2(self):
        B = -self.quiver2.matrix
        Qp = Quiver(B)
        for R in sink_set(self.quiver1):
            assert not any([Rp < Qp for Rp in isomorphismClass(R, permutations(3)) if (np.tril(Rp.matrix) >= 0).all()])

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

class MutationClassInitTests(unittest.TestCase):
    def setUp(self):
        #would be nice to support non-4 vertex quivers, especially small ones.
        self.classes = []
        self.acyclicL = []
        self.gcd_vecs = []
        self.members = []
        self.mat_ranks = []
        self.alexander_polys = []
        self.has_vortex = []
        Q = isolatedQuiver(4)
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(True)
        self.gcd_vecs.append((0,0,0,0))
        self.members.append(Q)
        self.mat_ranks.append(0)
        self.alexander_polys.append(polynomial.polynomial([1,-4,6,-4,1]))
        self.has_vortex.append(False)
        Q = Quiver(np.matrix([[0,-3,0,0],[3,0,0,0], [0,0,0,0], [0,0,0,0]]))
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(True)
        self.gcd_vecs.append((3,3,0,0))
        self.members.append(Q)
        self.mat_ranks.append(2)
        self.alexander_polys.append(polynomial.polynomial([1,-4,6,-4,1]) + 9*polynomial.polynomial([0,1,-2,1]))
        self.has_vortex.append(False)
        Q = Quiver(np.matrix([[0,3,2,0],[-3,0,3,0],[-2,-3,0,0], [0,0,0,0]]))
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(True)
        self.gcd_vecs.append((1,3,1,0))
        self.members.append(Q)
        self.mat_ranks.append(2)
        self.alexander_polys.append(polynomial.polynomial([-1,1])*polynomial.polynomial([-1,-37,37,1]))
        self.has_vortex.append(False)
        Q = Quiver(np.matrix([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]]))
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(False)
        self.gcd_vecs.append((2,2,2,0))
        self.members.append(Q)
        self.mat_ranks.append(2)
        self.alexander_polys.append(polynomial.polynomial([-1,1])*polynomial.polynomial([-1,-1,1,1]))
        self.has_vortex.append(Unknown) #in fact, True, if we use connectivity.
        R = Quiver(-np.matrix([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]]))
        self.classes.append(mutationClass(vertices=[R,Q], edges={R: {Q:1}, Q:{R:2}}, perms=permutations(4)))
        self.acyclicL.append(False)
        self.gcd_vecs.append((2,2,2,0))
        self.members.append(Q)
        self.mat_ranks.append(2)
        self.alexander_polys.append(polynomial.polynomial([-1,1])*polynomial.polynomial([-1,-1,1,1]))
        self.has_vortex.append(Unknown) #in fact, True, if we use connectivity.
        Q = Quiver(np.matrix([[0,3,-3,0],[-3,0,3,0],[3,-3,0,0], [0,0,0,0]]))
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(False)
        self.gcd_vecs.append((3,3,3,0))
        self.members.append(Q)
        self.mat_ranks.append(2)
        self.alexander_polys.append(polynomial.polynomial([-1,1])*polynomial.polynomial([-1,3,-3,1]))
        self.has_vortex.append(Unknown) #in fact, True, if we use connectivity.
        Q = Quiver(np.matrix([[0,2,-4,7],[-2,0,7,16],[4,-7,0,6], [-7,-16,-6,0]]))
        self.classes.append(mutationClass(Q, perms=permutations(4)))
        self.acyclicL.append(False)
        self.gcd_vecs.append((1,1,1,1))
        self.members.append(Q)
        self.mat_ranks.append(4)
        self.alexander_polys.append(False)
        self.has_vortex.append(True)

        
    def testMembership(self):
        for i in range(len(self.classes)):
            assert self.members[i] in self.classes[i], f"Missing quiver {self.members[i].matrix}"

    def testAcyclicInit(self):
        for i in range(len(self.classes)):
            assert self.classes[i].mutationAcyclic == self.acyclicL[i], f"Incorrect acyclicity {i}"

    def testGCDInit(self):
        for i in range(len(self.classes)):
            assert self.classes[i].gcdVector == self.gcd_vecs[i], f"Incorrect gcd_vector {i}"

    def testBRankInit(self):
        for i in range(len(self.classes)):
            assert self.classes[i].B_rank == self.mat_ranks[i], f"Incorrect matrix rank {i}"

    def testEq(self):
        assert self.classes[0] != self.classes[1]
        assert self.classes[1] != self.classes[2]
        assert self.classes[2] != self.classes[3]
        assert self.classes[2] != self.classes[4]
        assert self.classes[0] == self.classes[0]
        C = mutationClass(Quiver(np.matrix([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]])), perms=permutations(4))
        C.update()
        assert self.classes[3] == C, "Incorrect equality from update"
    
    def testAlexander(self):
        for i in range(len(self.classes)):
            assert self.classes[i].alexander_poly == self.alexander_polys[i], f"Incorrect alexander poly {i, str(self.classes[i].alexander_poly)}"
            assert self.classes[i].totallyProper == True if self.acyclicL[i] else Unknown

    def testHasVortexInit(self):
        for i in range(len(self.classes)):
            assert self.classes[i].hasVortex == self.has_vortex[i], f"incorrect vortex initialization for {i, self.classes[i].hasVortex}"

class updateMutationClassTests(unittest.TestCase):
    def setUp(self):
        self.markovClass = mutationClass(Quiver(np.matrix([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]])), perms=permutations(4))
        self.markovClass.update()
        self.markovResult = self.markovClass.update()

        self.isolatedClass = mutationClass(isolatedQuiver(4), perms=permutations(4))
        self.isolatedResult = self.isolatedClass.update()

        self.A2Class = mutationClass(Quiver(np.matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,0], [0,0,0,0]])), perms=permutations(4))
        self.A2Class.update()
        self.A2Result = self.A2Class.update()

        self.A3Class = mutationClass(Quiver(np.matrix([[0,1,0,0],[-1,0,1,0],[0,-1,0,0], [0,0,0,0]])), perms=permutations(4))
        self.A3Class.update()
        self.A3Class.update()
        self.A3Class.update()
        self.A3Class.update()
        self.A3Result = self.A3Class.update()

        self.AcyclicEventuallyClass = mutationClass(Quiver(np.matrix([[0,2,-3,0],[-2,0,2,0],[3,-2,0,0], [0,0,0,0]])), perms=permutations(4))
        self.AcyclicEventuallyClass.update()
        self.AcyclicEventuallyClass.update()
        self.AcyclicEventuallyClass.update()
        self.AcyclicEventuallyResult = self.AcyclicEventuallyClass.update()

        self.bigClass = mutationClass(Quiver(np.matrix([[0,2000,-300,0],[-2000,0,20,0],[300,-20,0,0], [0,0,0,0]])), perms=permutations(4))
        self.bigClass.update()
        self.bigClassResult = self.bigClass.update()

        self.vortexClass = mutationClass(Quiver(np.matrix([[0,2,-4,6],[-2,0,2,6],[4,-2,0,6], [-6,-6,-6,0]])), perms=permutations(4))
        self.vortexClass.update()
        self.vortexClass.update()
        self.vortexClassResult=self.vortexClass.update()

        self.vortexConnClass = mutationClass(Quiver(np.matrix([[0,2,-40,6],[-2,0,2,6],[40,-2,0,6], [-6,-6,-6,0]])), perms=permutations(4))
        self.vortexConnClass.update()
        self.vortexConnClass.update()
        self.vortexConnClassResult=self.vortexClass.update()
    
    def testResults(self):
        assert len(self.isolatedResult) == 0
        assert len(self.markovResult) == 0
        assert len(self.A3Result) == 0
        assert len(self.A2Result) == 0
        assert len(self.AcyclicEventuallyResult) != 0
        assert len(self.bigClassResult) == 0 #pruned away
        assert len(self.vortexClassResult) !=0


    def testUpdateMembership(self):
        assert Quiver(np.matrix([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]])) in self.markovClass
        assert isolatedQuiver(4) in self.isolatedClass
        assert isolatedQuiver(4) not in self.A3Class
        assert Quiver(np.matrix([[0,1,-1,0],[-1,0,1,0],[1,-1,0,0], [0,0,0,0]])) in self.A3Class
        assert Quiver(np.matrix([[0,-1,0,0],[1,0,1,0],[0,-1,0,0], [0,0,0,0]])) in self.A3Class
        assert Quiver(np.matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,0], [0,0,0,0]])) not in self.A3Class
        assert Quiver(np.matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,0], [0,0,0,0]])) in self.A2Class
        assert Quiver(np.matrix([[0,2,-3,0],[-2,0,2,0],[3,-2,0,0], [0,0,0,0]])).mutate(1).mutate(0).mutate(1).mutate(0) in self.AcyclicEventuallyClass
        assert Quiver(-np.matrix([[0,2,-5,0],[-2,0,6,0],[5,-6,0,0], [0,0,0,0]])) in self.AcyclicEventuallyClass
        assert Quiver(np.matrix([[0,2,-4,0],[-2,0,6,0],[4,-6,0,0], [0,0,0,0]])) not in self.AcyclicEventuallyClass
        assert Quiver(np.matrix([[0,2,-4,-6],[-2,0,2,-6],[4,-2,0,-6], [6,6,6,0]])) in self.vortexClass
        assert Quiver(np.matrix([[0,2,-4,-6],[-2,0,2,-6],[4,-2,0,-6], [6,6,6,0]])).mutate(1) in self.vortexClass
        

    def testUpdatedIntersection(self):
        assert self.markovClass.intersection(self.A2Class) is None
        assert self.markovClass.intersection(self.A3Class) is None
        assert self.markovClass.intersection(self.AcyclicEventuallyClass) is None
        assert self.markovClass.intersection(self.isolatedClass) is None
        assert self.markovClass.intersection(self.markovClass) == self.markovClass

    def testUpdatedUnion(self):
        assert self.markovClass.union(self.A2Class) is None
        assert self.markovClass.union(self.A3Class) is None
        assert self.markovClass.union(self.AcyclicEventuallyClass) is None
        assert self.markovClass.union(self.isolatedClass) is None
        assert self.markovClass.union(self.markovClass) == self.markovClass

    def testUpdatedFinite(self):
        assert self.markovClass.finite
        assert self.A2Class.finite
        assert self.isolatedClass.finite
        assert self.A3Class.finite
        assert not self.AcyclicEventuallyClass.finite
        assert not self.bigClass.finite
        assert not self.vortexClass.finite

    def testUpdatedConnected(self):
        assert not self.markovClass.mutationConnected 
        assert not self.A2Class.mutationConnected 
        assert not self.A3Class.mutationConnected 
        assert not self.isolatedClass.mutationConnected 
        assert self.vortexClass.mutationConnected is not False
        assert self.vortexConnClass.mutationConnected is not False

    def testLeastEdges(self):
        assert self.markovClass.leastEdges == 6
        assert self.A2Class.leastEdges == 1
        assert self.AcyclicEventuallyClass.leastEdges == 3
        assert self.A3Class.leastEdges == 2
        assert self.isolatedClass.leastEdges == 0

    def testRep(self):
        assert Quiver(np.matrix([[0,2,-3,0],[-2,0,2,0],[3,-2,0,0], [0,0,0,0]])) not in self.AcyclicEventuallyClass.possibleReps
        assert Quiver(np.matrix([[0,2,-1,0],[-2,0,0,0],[1,0,0,0], [0,0,0,0]])) in self.AcyclicEventuallyClass.possibleReps
        assert Quiver(np.matrix([[0,1,-1,0],[-1,0,1,0],[1,-1,0,0], [0,0,0,0]])) not in self.A3Class.possibleReps

    def testMuCyclic(self):
        assert not self.markovClass.mutationAcyclic
        assert self.markovClass.mutationCyclicSubquiver
        assert self.markovClass.mutationCyclicSubquiverWitness.hasMutCyclicSubquiver()
        assert self.A2Class.mutationAcyclic
        assert self.A3Class.mutationAcyclic
        assert not self.A2Class.mutationCyclicSubquiver
        assert self.A2Class.mutationCyclicSubquiverWitness is None
        assert self.AcyclicEventuallyClass.mutationAcyclic is Unknown
        assert not self.bigClass.mutationAcyclic

    def testForklessParts(self):
        assert self.markovClass.finitePFP
        assert self.markovClass.finite
        assert self.A2Class.finite
        assert self.A3Class.finiteFP
        assert self.bigClass.finiteFP is Unknown
        assert self.bigClass.finitePFP is Unknown

    def testBigWeightPruning(self):
        assert self.bigClass.hit_max_weight
        assert not self.A2Class.hit_max_weight
        assert not self.AcyclicEventuallyClass.hit_max_weight

    def testHasVortex(self):
        assert self.bigClass.hasVortex is Unknown, f"disconnected large fork incorrectly knows about vortices: {self.bigClass.hasVortex}"
        assert self.AcyclicEventuallyClass.hasVortex is False
        assert self.vortexClass.hasVortex
        assert self.vortexConnClass.hasVortex

    def testClassAbundance(self):
        assert not self.isolatedClass.mutationAbundant
        assert not self.bigClass.mutationAbundant
        assert not self.vortexClass.mutationAbundant
        assert self.vortexConnClass.mutationAbundant

    def testClassSurfaceQuiver(self):
        assert self.isolatedClass.is_surface_quiver
        assert self.A2Class.is_surface_quiver
        assert self.A3Class.is_surface_quiver
        assert not self.bigClass.is_surface_quiver
        assert not self.AcyclicEventuallyClass.is_surface_quiver

    def testAlexanderPreserved(self):
        for Q in self.A3Class.vertices:
            assert Q.alexander_poly() == self.A3Class.alexander_poly, f"alexander poly varies across class, {str(self.A3Class.alexander_poly)} vs {str(Q.alexander_poly())}"
        assert self.A3Class.totallyProper
        for Q in self.AcyclicEventuallyClass.vertices:
            assert Q.alexander_poly() == self.AcyclicEventuallyClass.alexander_poly, f"alexander poly varies across class, {str(self.AcyclicEventuallyClass.alexander_poly)} vs {str(Q.alexander_poly())}"
        assert self.AcyclicEventuallyClass.totallyProper
        assert not self.vortexClass.alexander_poly
        assert not self.vortexConnClass.alexander_poly
        assert not self.vortexClass.totallyProper
        assert not self.vortexConnClass.totallyProper

    def testCasals(self):
        c = self.A3Class.casals_det
        for Q in self.A3Class.vertices:
            assert c == Q.casals_det()
        c = self.vortexClass.casals_det
        for Q in self.vortexClass.vertices:
            assert c == Q.casals_det()

    def testSevenUnal(self):
        s = self.A2Class.sevens_ker
        for Q in self.A2Class.vertices:
            assert s == Q.seven_congruence()
        s = self.bigClass.sevens_ker
        for Q in self.bigClass.vertices:
            assert s == Q.seven_congruence()

        
class SurfaceQuiverTests(unittest.TestCase):
    def setUp(self):
        n = 6
        i,j,k,l,h,g = 0,1,2,3,4,5
        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,l] += 1
        self.BI = Quiver(B - np.transpose(B))
        self.BIm = self.BI.mutate(0)

        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[k,i] += 1
        self.BII = Quiver(B - np.transpose(B))
        self.BIIm = self.BII.mutate(0).mutate(1).mutate(2).mutate(1)

        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[j,i] += 1
        B[k,i] += 1
        self.BIIIa = Quiver(B - np.transpose(B))
        self.BIIIam = self.BIIIa.mutate(0)

        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,j] += 1
        B[i,k] += 1
        self.BIIIb = Quiver(B - np.transpose(B))

        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,j] += 1
        B[i,k] += 1
        B[j,l] += 1
        B[k,l] += 1
        B[l,i] += 1
        self.BIV = Quiver(B - np.transpose(B))

        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[j,h] += 1
        B[k,i] += 1
        B[h,i] += 1
        B[g,k] += 1
        B[g,h] += 1
        B[i,g] += 1
        self.BV = Quiver(B - np.transpose(B))
        self.BVm = self.BV.mutate(0)

        #non egs
        self.MI = Quiver([[0,3],[-3,0]])

        self.MII = Quiver([[0,2,0],[-2,0,2],[0,-2,0]])

        self.MIII = Quiver([[0,2,0],[-2,0,1],[0,-1,0]])

        self.MIV = Quiver([[0,0,1,1,1],[0,0,-1,-1,-1],[-1,1,0,0,0],[-1,1,0,0,0],[-1,1,0,0,0]])
        self.MIVm = self.MIV.mutate(0).mutate(2).mutate(3)

        #see SlowSurfaceTest...
        #self.E6 = Quiver([[0,1,0,0,0,0],[-1,0,1,0,0,0],[0,-1,0,-1,-1,0],[0,0,1,0,0,0],[0,0,1,0,0,1],[0,0,0,0,-1,0]])
        #self.X6 = Quiver([[0,1,-1,1,-1,-1],[-1,0,2,0,0,0],[1,-2,0,0,0,0],[-1,0,0,0,2,0],[1,0,0,-2,0,0],[1,0,0,0,0,0]])


    def testSurfaceBlocks(self):
        for Q in [self.BI, self.BII, self.BIIIa, self.BIIIb, self.BIV, self.BV]:
            r = surface_quiver(Q)
            assert r, f"couldn't recognize {Q.matrix}"
            blocks = r[1]
            count_zero = 0
            for b in blocks:
                assert all([0 <= x < Q.n for x in b[1]]), f"produced impossible index in block {Q.matrix, b}"
                if 0 in b[1]:
                    count_zero += 1
                assert len(b[1])==len(set(b[1])), f"block has repeated vertex {Q.matrix, b}"
            assert count_zero <= 2, f"index zero was too prevalent {Q.matrix, r}"
        for M in [self.MI, self.MII, self.MIII, self.MIV]:
            assert not surface_quiver(M), f"incorrectly found decomposition for {M.matrix} with decomposition {surface_quiver(M)[1]}"

    def testSurfaceBlocksMutated(self):
        for Q in [self.BIm, self.BIIm, self.BIIIam, self.BVm]:
            assert surface_quiver(Q), f"couldn't recognize {Q.matrix}"
        assert not surface_quiver(self.MIVm), f"incorrectly marked {self.MIVm.matrix} as being a surface quiver"

    def testSurfaceBlocksIso(self):
        p1 = [3,2,0,1,5,4]
        p2 = [5,4,3,2,1,0]
        for Q in [self.BI, self.BII, self.BIIIa, self.BV, self.BIIm, self.BVm]:
            assert surface_quiver(isomorphicQuiver(Q,p1)), f"couldn't recognize permuted block {p1, Q.matrix}"
            assert surface_quiver(isomorphicQuiver(Q,p2)), f"couldn't recognize permuted block {p2, Q.matrix}"

        for M in [self.MI, self.MII, self.MIII, self.MIV]:
            assert not surface_quiver(isomorphicQuiver(M,[1,0,2,3,4,5])), f"incorrectly recognizes permuted block {[1,0,2,3,4,5], M.matrix} as {surface_quiver(isomorphicQuiver(M,p1))}"
            assert not surface_quiver(isomorphicQuiver(M,[1,0,2,4,3,5])), f"incorrectly recognizes permuted block {[1,0,2,4,3,5], M.matrix}"

@unittest.skip("Larger surface quiver tests skipped.")
class SlowSurfaceTests(unittest.TestCase):
    def setUp(self):
        n = 6
        i,j,k,l,h,g = 0,1,2,3,4,5
        B = np.matrix(np.zeros((n,n), dtype='object'), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[j,h] += 1
        B[k,i] += 1
        B[h,i] += 1
        B[g,k] += 1
        B[g,h] += 1
        B[i,g] += 1
        self.BV = Quiver(B - np.transpose(B))
        self.BVm = self.BV.mutate(0)

        i,j,k,l,h,g = 0,1,2,3,4,5
        B = np.matrix(np.zeros((8,8), dtype='object'), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[j,h] += 1
        B[k,i] += 1
        B[h,i] += 1
        B[g,k] += 1
        B[g,h] += 1
        B[i,g] += 1
        j,k,l = 3,6,7
        B[i,j] += 1
        B[i,k] += 1
        B[j,l] += 1
        B[k,l] += 1
        B[l,i] += 1
        self.BVp = Quiver(B - np.transpose(B))
        self.BVpm = self.BV.mutate(0)

        #non egs
        self.E6 = Quiver([[0,1,0,0,0,0],[-1,0,1,0,0,0],[0,-1,0,-1,-1,0],[0,0,1,0,0,0],[0,0,1,0,0,1],[0,0,0,0,-1,0]])
        self.X6 = Quiver([[0,1,-1,1,-1,-1],[-1,0,2,0,0,0],[1,-2,0,0,0,0],[-1,0,0,0,2,0],[1,0,0,-2,0,0],[1,0,0,0,0,0]])


    def testSurfaceBlocks(self):
        for Q in [self.BV, self.BVp, self.BVpm]:
            r = surface_quiver(Q)
            assert r, f"couldn't recognize {Q.matrix}"
            blocks = r[1]
            count_zero = 0
            for b in blocks:
                assert all([0 <= x < Q.n for x in b[1]]), f"produced impossible index in block {Q.matrix, b}"
                if 0 in b[1]:
                    count_zero += 1
                assert len(b[1])==len(set(b[1])), f"block has repeated vertex {Q.matrix, b}"
            assert count_zero <= 2, f"index zero was too prevalent {Q.matrix, r}"
        for M in [self.E6, self.X6]:
            assert not surface_quiver(M), f"incorrectly found decomposition for {M.matrix} with decomposition {surface_quiver(M)[1]}"

    def testSurfaceBlocksMutated(self):
        for Q in [self.BVm]:
            assert surface_quiver(Q), f"couldn't recognize {Q.matrix}"
        assert not surface_quiver(self.E6.mutate(1)), f"incorrectly marked {self.E6.matrix} as being a surface quiver"

    def testSurfaceBlocksIso(self):
        p1 = [3,2,0,1,5,4,6,7]
        p2 = [5,4,3,2,1,0,6,7]
        for Q in [self.BV, self.BVp, self.BVpm]:
            assert surface_quiver(isomorphicQuiver(Q,p1)), f"couldn't recognize permuted block {p1, Q.matrix}"
            assert surface_quiver(isomorphicQuiver(Q,p2)), f"couldn't recognize permuted block {p2, Q.matrix}"

        for M in [self.E6, self.X6]:
            assert not surface_quiver(isomorphicQuiver(M,[1,0,2,3,4,5])), f"incorrectly recognizes permuted block {[1,0,2,3,4,5], M.matrix} as {surface_quiver(isomorphicQuiver(M,p1))}"
            assert not surface_quiver(isomorphicQuiver(M,[1,0,2,4,3,5])), f"incorrectly recognizes permuted block {[1,0,2,4,3,5], M.matrix}"


class TestLinAlg(unittest.TestCase):
    def setUp(self):
        self.quiver = Quiver(np.matrix([[0,3,-3,0],[-3,0,3,0],[3,-3,0,0], [0,0,0,0]]))
        self.M = np.matrix([[2,0],[0,3]])
        self.N = np.matrix([[2,-1],[-2,2]])
        self.vortex = Quiver(np.matrix([[0,1,2,3],[-1,0,2,-2],[-2,-2,0,3],[-3,2,3,0]]))

    def testEign(self):
        np.testing.assert_almost_equal(eigenvalues(self.quiver.matrix), np.array([0+3*np.emath.sqrt(-3.), 0-3*np.emath.sqrt(-3.), 0, 0]))
        np.testing.assert_almost_equal(eigenvalues(self.M), np.array([2.,3.]))
        np.testing.assert_almost_equal(eigenvalues(self.N), np.array([2.+np.sqrt(2.), 2.-np.sqrt(2.)]))


    def testRank(self):
        assert 2 == matrix_rank(self.quiver.matrix), f"incorrect rank of disconnected markov"
        assert 2 == matrix_rank(self.M)
        assert 0 == matrix_rank(isolatedQuiver(4).matrix)
        assert 4 == matrix_rank(self.vortex.matrix)

class IsomorphismClassTests(unittest.TestCase):
    def setUp(self):
        self.markov = Quiver(np.matrix([[0,2,-2],[-2,0,2],[2,-2,0]]))
        self.empty = Quiver(np.matrix([[0]]))
        self.A3 = Quiver(np.matrix([[0,1,0],[-1,0,1],[0,-1,0]]))

    def testNotIso(self):
        assert self.markov not in isomorphismClass(self.empty, permutations(1))
        assert self.A3 not in isomorphismClass(self.markov, permutations(3))
        assert self.A3.mutate(1) not in isomorphismClass(self.A3, permutations(3))
        assert self.A3 not in isomorphismClass(self.A3.mutate(1), permutations(3))
        assert self.A3.mutate(0) not in isomorphismClass(self.A3, permutations(3))
        assert self.markov not in isomorphismClass(self.A3, permutations(3))

    def testIso(self):
        assert self.markov in isomorphismClass(self.markov, permutations(3))
        assert self.empty in isomorphismClass(self.empty, permutations(1))
        assert self.markov.mutate(1) in isomorphismClass(self.markov, permutations(3))
        assert self.A3.mutate(0).mutate(2) in isomorphismClass(self.A3, permutations(3))
        assert self.A3.mutate(1).mutate(2) in isomorphismClass(self.A3, permutations(3))

class MutationCycleTests(unittest.TestCase):
    def setUp(self):
        self.A3Class = mutationClass(Quiver(np.matrix([[0,1,0,0],[-1,0,1,0],[0,-1,0,0], [0,0,0,0]])), perms=permutations(4))
        self.A3Class.update()
        self.A3Class.update()
        self.A3Class.update()
        self.A3Class.update()
        self.A3Class.update()

        self.markovClass = mutationClass(Quiver([[0,2,-2,0],[-2,0,2,0],[2,-2,0,0], [0,0,0,0]]), perms=permutations(4))
        self.markovClass.update()
        self.markovClass.update()
        self.markovClass.update()

        self.finiteForklessClass = mutationClass(Quiver([[0,19,19,-19],[-19,0,19,19],[-19,-19,0,19],[19,-19,-19,0]]), perms=permutations(4))
        self.finiteForklessClass.update()
        self.finiteForklessClass.update()

        self.A1Class = mutationClass(isolatedQuiver(1), perms=permutations(1))
        self.A1Class.update()

    def testMutCycleDetected(self):
        assert self.A1Class.hasMutationCycle, f"mutation cycle of length 1 was ignored"
        assert self.A3Class.hasMutationCycle
        assert self.markovClass.hasMutationCycle
        assert not self.finiteForklessClass.hasMutationCycle

    def testFiniteForklessDetected(self):
        assert self.A1Class.finite
        assert self.A1Class.finiteFP
        assert self.A1Class.finitePFP
        assert self.A3Class.finite
        assert self.A3Class.finiteFP
        assert self.A3Class.finitePFP
        assert self.markovClass.finite
        assert self.markovClass.finitePFP
        assert self.markovClass.finiteFP
        assert self.finiteForklessClass.finiteFP
        assert self.finiteForklessClass.finitePFP
        assert not self.finiteForklessClass.finite

class UnknownTests(unittest.TestCase):
    def testIs(self):
        assert Unknown is Unknown
        assert not Unknown is False
        assert not Unknown is True
        assert not False is Unknown
        assert not True is Unknown

    def testIsNot(self):
        assert Unknown is not False
        assert Unknown is not True
        assert False is not Unknown
        assert True is not Unknown
        assert not Unknown is not Unknown

    def testAnd(self):
        assert Unknown & True is Unknown
        assert Unknown & False is False
        assert Unknown & Unknown is Unknown
        assert (True & Unknown) is (Unknown & True)
        assert (False & Unknown) is (Unknown & False)

    def testOr(self):
        assert Unknown | True is True
        assert Unknown | False is Unknown
        assert Unknown | Unknown is Unknown
        assert (True | Unknown) is (Unknown | True)
        assert (False | Unknown) is (Unknown | False)

if __name__ == "__main__":
    unittest.main()
