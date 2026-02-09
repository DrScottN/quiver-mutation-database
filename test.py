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
        
class QuiverInvariants3TestCase(unittest.TestCase):
    def setUp(self):
        self.quiver_c1 = Quiver([[0,-3,6,0], [3,0,-5,0], [-6,5,0,0], [0,0,0,0]])
        self.quiver_c2 = Quiver([[0,3,-6,0], [-3,0,5,0], [6,-5,0,0], [0,0,0,0]])
        self.quiver_a2 = Quiver([[0,2,-6,0], [-2,0,5,0], [6,-5,0,0], [0,0,0,0]])
        self.quiver_a1 = Quiver([[0,-2,6,0], [2,0,-5,0], [-6,5,0,0], [0,0,0,0]])

    def testAcyclic(self):
        assert not self.quiver_a1.acyclic()
        assert not self.quiver_a2.acyclic()
        assert not self.quiver_c1.acyclic()
        assert not self.quiver_c2.acyclic()
    
    def testConnected(self):
        assert not self.quiver_a1.connected()
        assert not self.quiver_a2.connected()
        assert not self.quiver_c1.connected()
        assert not self.quiver_c2.connected()

    def testCyclicSubquiver(self):
        assert not self.quiver_a1.hasMutCyclicSubquiver(), 'mutation acyclic quiver is incorrectly marked as mutation cyclic'
        assert not self.quiver_a2.hasMutCyclicSubquiver(), 'mutation acyclic quiver is incorrectly marked as mutation cyclic'
        assert self.quiver_c1.hasMutCyclicSubquiver(), 'mutation cyclic subquiver is incorrectly unnoticed'
        assert self.quiver_c2.hasMutCyclicSubquiver(), 'mutation cyclic subquiver is incorrectly unnoticed'

    def testCycles(self):
        assert [0,1,2] in list(self.quiver_a1.chordless_cycles()) or [0,2,1] in list(self.quiver_a1.chordless_cycles()), 'missing chordless cycle'
        assert [0,1,2] in list(self.quiver_a2.chordless_cycles()) or [0,2,1] in list(self.quiver_a2.chordless_cycles()), 'missing chordless cycle'
        assert [0,1,2] in list(self.quiver_c1.chordless_cycles()) or [0,2,1] in list(self.quiver_c1.chordless_cycles()), 'missing chordless cycle'
        assert [0,1,2] in list(self.quiver_c2.chordless_cycles()) or [0,2,1] in list(self.quiver_c2.chordless_cycles()), 'missing chordless cycle'
        assert len(list(self.quiver_a1.chordless_cycles()))==1, 'incorrectly found more cycles'
        assert len(list(self.quiver_a2.chordless_cycles()))==1, 'incorrectly found more cycles'
        assert len(list(self.quiver_c1.chordless_cycles()))==1, 'incorrectly found more cycles'
        assert len(list(self.quiver_c2.chordless_cycles()))==1, 'incorrectly found more cycles'


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
        
    def testSetChanges(self):
        assert self.quiver.updateWeight(0,0,1).matrix[0][1] == 0, "update weight does not update the weight"
        assert hash(self.quiver.updateWeight(0,0,1)) != hash(self.quiver_dup), f"setting a weight doesn't change hash"
        assert hash(self.quiver.updateWeight(0,0,1)) != hash(self.quiver), f"setting a weight doesn't change hash"


if __name__ == "__main__":
    unittest.main()