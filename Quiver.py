import copy
from collections import Counter
import numpy as np 
import itertools
import math

import polynomial

class Quiver():
    def __init__(self, matrix, validate = True):
        if validate and not isinstance(matrix, np.matrix):
            if len(matrix) != len(matrix[0]):
                raise Exception("Not a square matrix")
            else:
                for i, l in enumerate(matrix):
                    for j, a in enumerate(l):
                        if a != -matrix[j][i]:
                            raise Exception("Not a skew-symmetric matrix")

        self.matrix = np.matrix(matrix, dtype=object) #matrix # type object lets us put arbitrary python in there.
        self.n = len(matrix)
        self.vertices = [v for v in range(self.n)]
        self.numEdges = sum(self.matrix[i,j] for i in self.vertices for j in self.vertices if self.matrix[i,j] > 0)

    def __hash__(self):
        t = tuple(self.matrix[i,j] for i in range(self.n) for j in range(i+1,self.n))
        return hash(t)

    def __str__(self):
        return '\n'.join(str(r) for r in self.matrix)

    def __eq__(self, other):
        return (self.matrix == other.matrix).all() #previously no .all

    def __lt__(self, other):
        # Lexicographical order on matrix elements
        if self.__eq__(other):
            return False
        elif self.n != other.n:
            raise Exception("Attempted to compare two quivers of different size")

        for i in range(self.n):
            for j in range(self.n):
                if self.matrix[i,j] < other.matrix[i,j]:
                    return True
                elif self.matrix[i,j] > other.matrix[i,j]:
                    return False

        return False

    def __deepcopy__(self, memo):
        return Quiver(copy.deepcopy(self.matrix), validate=False)

    def determinant(self):
        if self.n % 2 == 1:
            return 0
        
        def pfaffian(matrixAsList):
            n = len(matrixAsList)

            if n == 0:
                return 1
            
            #print(n)
            return sum([pow(-1,j+1) * matrixAsList[0][j] * pfaffian([[matrixAsList[i][k] for i in range(1,n) if i != j] for k in range(1,n) if k != j]) for j in range(1,n)])

        return pfaffian(self.matrix.tolist())**2
    
    def subquiverRemoveOneVertex(self, v):
        # Takes the subquiver by removing the vertex v
        return Quiver([[self.matrix[i,j] for i in range(self.n) if i != v] for j in range(self.n) if j != v])

    def connected(self):
        # Finds if the quiver is connected as a simple graph
        seen = [0] + [j for j in range(1,self.n) if self.matrix[0,j] != 0]
        numSeen = 1

        while numSeen < len(seen):
            extension = []
            for j in seen[numSeen:]:
                extension.extend([k for k in range(self.n) if self.matrix[j,k] != 0 and k not in extension])

            numSeen = len(seen)
            seen.extend([k for k in extension if k not in seen])

        return numSeen == self.n
    
    def sources(self):
        # Get the sources in a quiver

        return [i for i in range(self.n) if all([self.matrix[i,j]>=0 for j in range(self.n)])]
    
    def sinks(self):
        # Get the sources in a quiver

        return [i for i in range(self.n) if all([self.matrix[i,j]<=0 for j in range(self.n)])]

    def acyclic(self):
        # Finds if the quiver is acyclic or not
        if self.n < 3:
            return True

        #sources = [i for i in range(self.n) if sum(self.matrix[i]) == sum([abs(self.matrix[i,j]) for j in range(self.n)])]
        sources = self.sources() 

        if len(sources) == 0:
            return False
        elif len(sources) == self.n:
            return True
        
        vertices = [i for i in range(self.n) if i not in sources]

        newMatrix = [[self.matrix[i,j] for i in vertices] for j in vertices]

        return Quiver(newMatrix).acyclic()

    def complete(self):
        # Finds if the quiver is complete or not
        terms = [self.matrix[i,j] for i in range(self.n) for j in range(i+1,self.n)]

        s = 1

        for t in terms:
            s *= t

        return s != 0

    def abundant(self):
        # Finds if the quiver is abundant
        for i in range(self.n):
            for j in range(i+1, self.n):
                if abs(self.matrix[i,j]) < 2:
                    return False

        return True 
    
    def forkWithPOR(self, r):
        # returns true if a fork with point of return r
        # else false
        if self.acyclic():
            return False
        elif not self.abundant(): 
            return False
        elif not self.subquiverRemoveOneVertex(r).acyclic():
            return False
        
        for i in range(self.n):
            if i == r:
                continue

            if self.matrix[i,r] > 0:
                for j in range(self.n):
                    if j == r or j == i:
                        continue

                    if self.matrix[r,j] > 0:
                        value = self.matrix[j,i] > 0
                        value = value and self.matrix[i,r] < self.matrix[j,i]
                        value = value and self.matrix[r,j] < self.matrix[j,i]

                        if not value:
                            return False
                        
            if self.matrix[i,r] < 0:
                for j in range(self.n):
                    if j == r or j == i:
                        continue

                    if self.matrix[r,j] < 0:
                        value = self.matrix[i,r] > self.matrix[j,i]
                        value = value and self.matrix[r,j] > self.matrix[j,i]
                        value = value and self.matrix[j,i] < 0

                        if not value:
                            return False

        return True
   
    def preForkWithVertices(self, r, i, j):
        # returns true if pre-fork with given vertices
        # else false
        if self.acyclic():
            return False
        
        Q = self.subquiverRemoveOneVertex(i)
        P = self.subquiverRemoveOneVertex(j)

        if not Q.forkWithPOR(r) or not P.forkWithPOR(r):
            return False

        for k in self.vertices:
            if self.matrix[i,k]*self.matrix[k,j] > 0:
                return False

        return True

    def preForkWithPOR(self, r):
        # Returns true if pre-fork with point of return r
        # else false (even if it is a fork with por r)
        
        for i in self.vertices:
            for j in self.vertices:
                if j == i or r in [i,j]:
                    continue

                if self.matrix[i,j] > 1:
                    continue

                if self.preForkWithVertices(r, i, j):
                    return True

        return False

    def threeCycle(self):
        # Returnes whether the quiver is a 3-cycle or not
        if self.n != 3:
            return False
        
        return not self.acyclic()

    def hasMutCyclicSubquiver(self):
        # Returns whether the quiver has a 3-cycle as a subquiver
        if self.n != 4:
            raise Exception("Not implemented for more than four vertices")
        
        def markovInvariant(matrix):
            if len(matrix) != 3:
                raise Exception("Not implemented yet")
            
            return (matrix[0,1]**2) + (matrix[1,2]**2) + (matrix[2,0]**2) - abs(matrix[0,1]*matrix[1,2]*matrix[2,0])

        return any(Q.threeCycle() and markovInvariant(Q.matrix) <= 4 and Q.abundant() for Q in [self.subquiverRemoveOneVertex(i) for i in self.vertices])

    def vortex(self):
        # Returns whether the quiver is a vortex or not
        if self.n != 4:
            return False
        
        vertices = self.sources() + self.sinks()

        if len(vertices) != 1:
            return False
        
        vertex = vertices[0]
        
        if 0 in self.matrix[vertex,:vertex] or 0 in self.matrix[vertex,vertex+1:]:
            return False
        
        Q = self.subquiverRemoveOneVertex(vertex)

        return Q.threeCycle()

    def oppositeQuiver(self):
        #return Quiver([[-self.matrix[j,i] for i in range(self.n)] for j in range(self.n)])
        return Quiver(-1*self.matrix)
    
    def mutate(self, k):
        # Gives the quiver formed by mutating at vertex k

        newMatrix = copy.deepcopy(self.matrix)

        for i in range(self.n):
            for j in range(self.n):
                if k in [i,j]:
                    newMatrix[i,j] = -newMatrix[i,j]
                elif self.matrix[i,k]*self.matrix[k,j] > 0:
                    add = self.matrix[i,k]*self.matrix[k,j] if self.matrix[i,k] > 0 else -self.matrix[i,k]*self.matrix[k,j]
                    newMatrix[i,j] += add

        return Quiver(newMatrix)

    def hasSourceSink(self):
        # Checks whether a source or sink exists

        return len(self.sources()) + len(self.sinks()) > 0     

    def chordless_cycles(self):
        # Iterator for all chordless cycles in the quiver. Formatted as list of vertex indices.
        
        def neighbors(v):
            return [j for j in range(self.n) if self.matrix[v,j] != 0]
        
        def chordless_cycles_at(v):
            # Iterator for chordless cycles with minimal element v
            paths = []
            ends = [u for u in neighbors(v) if u > v]
            for u in ends[:-1]:
                paths.append([[v,u], [v,u]])
            
            while len(paths)!=0:
                P, chord = paths.pop()
                N = [u for u in neighbors(P[-1]) if u not in chord and u > v]
                for u in N:
                    if u in ends:
                        if u > P[1]: #deduplication
                            yield P+[u]
                        continue
                    paths.append([P+[u], chord + N])   

        for i in range(self.n):
            for C in chordless_cycles_at(i):
                yield C

    def winding_data(self, sigma):
        # Using linear order sigma (a permutation of the indices), return winding numbers and edges against direction of traversal of the chordless cycles.
        winding_dict = dict()
        for C in self.chordless_cycles():
            wind = 0 #times orbited
            lefts = 0 #arrows against flow
            for i in range(len(C)):
                wind += (sigma[C[(i+1)%len(C)]] < sigma[C[i]])
                lefts += self.matrix[C[i],C[(i+1)%len(C)]] < 0
            winding_dict[tuple(C)] = (wind - lefts, lefts)
        return winding_dict


    def cyclic_order(self):
        # Return the cyclic ordering of the quiver with all chordless cycles having winding number 1 (if oriented) or 0 (if acyclic). 
        #  Returns False if no such order exists.
        if self.n==1:
            return [0]
        if self.n==2:
            return [0,1]
        for sigma_p in permutations(self.n-1): #[[0,1,2,3], [0,2,1,3], [0,1,3,2], [0,2,3,1], [0,3,1,2], [0,3,2,1]]:
            sigma = sigma_p + [self.n-1]
            W = self.winding_data(sigma)
            good = True
            for C in W.keys():
                wind, left = W[C]
                if left > len(C)/2:
                    wind = -wind
                    left = len(C) - left
                if wind not in [0,1]:
                    good=False
                    break
                if wind==0 and left ==0:
                    good=False
                    break
                if wind==1 and left > 0:
                    good=False
                    break
            if not good:
                continue
            return sigma
        return False

    def Umatrix(self):
        # Compute a unipotent companion. 
        #  Return False if this quiver has no potentailly totally proper order.
        
        sigma = self.cyclic_order()
        if sigma==False:
            return False
        U = np.matrix([[-(i<j)*self.matrix[sigma[i],sigma[j]] for j in range(self.n)] for i in range(self.n)], dtype=object)
        for i in range(self.n):
            U[i,i]=1
        return U


    def markov(self):
        #computes the markov invariant. Return False if this quiver has no potentially totally proper order.
        U = self.Umatrix()
        if U is False:
            return False
        #compute trace of U U^-1 (= U (I + N + N^2 + ...))
        N = np.transpose(np.matlib.identity(self.n, dtype=object) - U)
        Cosquare = U.copy()
        for e in range(1,self.n):
            Cosquare += U * N**e
        trace = 0
        for i in range(self.n):
            trace += Cosquare[i,i]
        return self.n - trace
        
    def alexander_poly(self):
        # Computes the alexander polynomial (expressed as a polynomial object). Returns False if the quiver is not potentially totally proper.
        U = self.Umatrix()
        if U is False:
            return False
        Ux = np.matrix([[U[i,j]*polynomial.polynomial([0,1]) for i in range(self.n)] for j in range(self.n)], dtype=object) - np.matrix([[U[j,i]*polynomial.polynomial([1]) for i in range(self.n)] for j in range(self.n)], dtype=object)
        det = polynomial.polynomial([0])
        for p in permutations(self.n):
            a = polynomial.polynomial([1])
            for i in range(self.n):
                a *= Ux[i, p[i]]
            det += a * sign_perm(p)
        return det

    def gcd_vector(self):
        # Computes the gcd vector of self, as a tuple with ith entry the gcd of row i.
        return tuple(math.gcd(*[self.matrix[i,j] for j in range(self.n)]) for i in range(self.n))
        


def eigenvalues(M):
    # Computes the eigenvalues of the given matrix, cast to float.
    #   Note that the rank is a mutation invariant but not these values.
    return np.linalg.eig(np.matrix(self.matrix))


def sink_set(quiver):
    # Calculates all sink/source mutation equivalent quivers
    todo = [quiver]
    visited = []
    while todo:
        q = todo.pop()
        visited.append(q)
        for v in q.sinks():
            p = q.mutate(v)
            if p in visited:
                continue
            todo.append(p)
    return visited

def generate_acyclics_from_Alexander(alex):
    # Takes an Alexander polynomial alex and generates all acyclic quivers with that polynomial.
    #  Iterator. Assumes alex is uni-variate. Not efficient.
    degree = max([j[0] for j in alex.coeffDict.keys() if alex.coeffDict[j]!=0])
    markov = degree + alex.coeffDict[(degree-1,)]
    if markov < 0:
        return

    for w in itertools.product(list(range(math.isqrt(markov)+1)), repeat=(degree*(degree-1))//2):
        # iterates over all combinations of weights < sqrt(markov), which is much too large
        #  faster is to build subquivers that are smaller than markov. Even faster potentially is to use more coefficients and a solver.
        
        X = np.matrix([[0 for i in range(degree)] for j in range(degree)])
        #print(degree, np.tril_indices(X.shape[0], k = -1), w)
        X[np.tril_indices(X.shape[0], k = -1)] = w
        Q = Quiver(X - np.transpose(X))
        if Q.alexander_poly() != alex:
            continue
        #to avoid repeats; this is very slow.
        if any([any([R < Q for R in isomorphismClass(Rp, permutations(degree)) if (np.tril(R.matrix) >= 0).all()]) for Rp in sink_set(Q)]):
            continue

        yield Q


class mutationClass():
    # Object to hold mutation-equivalent quivers
    # Each quiver is represented as a vertex of the mutation class
    # Each mutation is represented as an edge of the mutation class
    # Edges are given as a dict of dicts. 
    # For example, for a quiver Q one mutation away from P at vertex k:
    #   {Q : {P : k}, P : {Q : k}}
    # Could easily use this to build a graph visualization of mutation classes

    # The class will be split up into two main ideas: essentially a set for mutations and a tool for exploring mutation-classes

    # Will implement mutationClass unions and intersections

    def __init__(self, Q = None, perms = None, vertices = None, edges = None):
        # Takes in a quiver Q as the first member of our mutation class. perms is required.
        if Q is None:
            if vertices is None or edges is None:
                raise Exception("No vertices or edges given. Try again")

            self.initializationFromVerticesAndEdges(vertices, edges)
            self.vertices.sort()
            self.initialQ = self.vertices[0]
        else:
            self.initialQ = copy.deepcopy(Q) 
            self.vertices = [self.initialQ]
            self.forefront = [self.initialQ] # These are all the quivers that have not been fully mutated
            self.edges = {self.initialQ : dict()} # self.edges[Q][P] gives the label of the mutation taking one to the other
            self.numVertices = self.initialQ.n
        
        self.perms = perms
        self.possibleReps = [self.initialQ]
        self.possibleIsoReps = [isomorphismClass(self.initialQ, perms)[0]]

        # Now, check for the mutation-invariants we can immediately find
        self.mutationConnected = self.initialQ.connected()
        self.determinant = self.initialQ.determinant()
        self.leastEdges = self.initialQ.numEdges
        self.size = 1
        self.gcdVector = self.initialQ.gcd_vector()


        if self.initialQ.acyclic():
            self.mutationAcyclic = True
            self.mutationCyclicSubquiver = False
            self.mutationCyclicSubquiverWitness = None
            self.couldBeMutationCyclic = False
        else:
            if self.initialQ.hasMutCyclicSubquiver():
                self.foundCyclicSubquiver(self.initialQ)
            else:
                self.mutationCyclicSubquiver = None
                self.mutationCyclicSubquiverWitness = None
                self.mutationAcyclic = None
            self.couldBeMutationCyclic = True

        if not self.initialQ.complete():
            self.mutationComplete = False
            self.couldBeMutationComplete = False
        else:
            self.mutationComplete = None
            self.couldBeMutationComplete = True

        if self.initialQ.vortex(): # need to update this for quivers without 4 vertices
            self.hasVortex = True
            self.couldBeMutationVortexFree = False
        else:
            self.hasVortex = None
            self.couldBeMutationVortexFree = True

        # Setup possible mutation-invariants
        self.finite = None
        self.finiteFP = None
        self.finitePFP = None

        self.couldBeFinite = True
        self.couldBeFiniteFP = True
        self.couldBeFinitePFP = True

    def __str__(self):
        return str(self.representative())

    
    def __hash__(self):
        return hash(self.representative())
    # Need a better hash method for this class

    def __eq__(self, other):
        return not self.emptyIntersection(other)

    def __contains__(self, Q):
        if not isinstance(Q, Quiver):
            return False

        return Q in self.vertices

    def __iter__(self):
        return iter(self.vertices)

    def initializationFromVerticesAndEdges(self, vertices, edges):
        # Allows us to initialize a mutationClass without starting with a single quiver
        # We assume that vertices lists all vertices
        # Edges may not have the edges for both sides
        # Edges does not include vertices not present in vertices

        self.vertices = copy.deepcopy(vertices) # Why the deepcopy? Want to be absolutely sure I'm not messing with unexpected things
        self.edges = copy.deepcopy(edges)
        self.forefront = []
        self.numVertices = self.vertices[0].n
        
        # Fix any messed up edges
        for Q in self.edges:
            for P in self.edges[Q]:
                if P not in self.edges:
                    self.edges[P] = {Q : self.edges[Q][P]}
                elif Q not in self.edges[P]:
                    self.edges[P][Q] = self.edges[Q][P]

        # Now can populate the forefront
        for Q in self.vertices:
            if Q not in self.edges:
                self.edges[Q] = dict()
            
            if len(self.edges[Q]) < self.numVertices:
                self.forefront.append(Q)

    def union(self, *others):
        """This returns the union of two mutation classes
            ONLY IF THEY HAVE A NONTRIVIAL INTERSECTION"""
        other = others[0]

        if len(others) > 1:
            other = other.union(*others[1:])

        if other is None:
            return None # Why? Want the union of two non-intersecting mutation-classes to be empty
        elif isinstance(other, Quiver):
            other = mutationClass(other, self.perms)

        # Now we are looking at the union of exactly two mutation classes

        selfVertices = set(self.vertices)
        otherVertices = set(other.vertices)
        vertices = list(selfVertices | otherVertices)
        commonVertices  = list(selfVertices & otherVertices)
        
        if len(commonVertices) == 0:
            # Again, this is returning an empty union if they do not intersect anywhere
            return None
        
        # Can now assume that they intersect at least one vertex

        selfEdges = {Q : {P : k for P, k in self.edges[Q].items() if P in vertices} for Q in commonVertices}
        otherEdges = {Q : {P : k for P, k in other.edges[Q].items() if P in vertices} for Q in commonVertices}
        commonEdges = {Q : {**selfEdges[Q], **otherEdges[Q]} for Q in commonVertices}
        commonEdges.update({Q : self.edges[Q] for Q in selfVertices if Q not in commonVertices})
        commonEdges.update({Q : other.edges[Q] for Q in otherVertices if Q not in commonVertices})

        newMutClass = mutationClass(vertices = vertices, edges = commonEdges, perms = self.perms)
        newMutClass.couldBeMutationCyclic = self.couldBeMutationCyclic and other.couldBeMutationCyclic
        newMutClass.couldBeMutationComplete = self.couldBeMutationComplete and other.couldBeMutationComplete
        newMutClass.couldBeMutationVortexFree = self.couldBeMutationVortexFree and other.couldBeMutationVortexFree
        newMutClass.couldBeFinite = self.couldBeFinite and other.couldBeFinite
        newMutClass.couldBeFiniteFP = self.couldBeFiniteFP and other.couldBeFiniteFP
        newMutClass.couldBeFinitePFP = self.couldBeFinitePFP and other.couldBeFinitePFP

        return newMutClass

    def intersection(self, *others):
        # Gets the intersection of two mutation classes
        # Returns a mutation class containing the common quivers
        # Else returns none if no quivers are common
        other = others[0]
        if len(others) > 1:
            other = other.intersection(*others[1:])

        if other is None:
            return None
        elif isinstance(other, Quiver):
            other = mutationClass(other, self.perms)

        # Now just have one intersection to do with the mutation class other

        selfVertices = set(self.vertices)
        otherVertices = set(other.vertices)
        verticesSet = selfVertices & otherVertices
        vertices  = list(verticesSet)

        if len(vertices) == 0:
            return None

        # Have a non-trivial intersection now

        selfEdges = {Q : {P : k for P, k in self.edges[Q].items() if P in vertices} for Q in vertices}
        otherEdges = {Q : {P : k for P, k in other.edges[Q].items() if P in vertices} for Q in vertices}
        edges = {Q : {**selfEdges[Q], **otherEdges[Q]} for Q in vertices}

        newMutClass = mutationClass(vertices = vertices, edges = edges, perms = self.perms)
        newMutClass.couldBeMutationCyclic = self.couldBeMutationCyclic and other.couldBeMutationCyclic
        newMutClass.couldBeMutationComplete = self.couldBeMutationComplete and other.couldBeMutationComplete
        newMutClass.couldBeMutationVortexFree = self.couldBeMutationVortexFree and other.couldBeMutationVortexFree
        newMutClass.couldBeFinite = self.couldBeFinite and other.couldBeFinite
        newMutClass.couldBeFiniteFP = self.couldBeFiniteFP and other.couldBeFiniteFP
        newMutClass.couldBeFinitePFP = self.couldBeFinitePFP and other.couldBeFinitePFP

        return newMutClass

    def emptyIntersection(self,other):
        return self.intersection(other) is None

    def representative(self):
        # Gets a quiver with lowest edge weights to represent the class
        if len(self.possibleReps) > 1:
            self.possibleReps.sort()
        
        return self.possibleReps[0]

    def isomorphicRepresentative(self):
        # Gets an isomorphic quiver with lowest edge weights to represent the class
        if len(self.possibleIsoReps) > 1:
            self.possibleIsoReps.sort()

        return self.possibleIsoReps[0]

    def update(self):
        # updates the mutation class by mutating every quiver on the edge in every possible direction
        # Naturally throws away forks and pre-forks if they have the same point of return as where we mutated
        newForefront = []

        for Q in self.forefront:
            for k in Q.vertices:
                P = Q.mutate(k) # Mutate the Quiver
                if P.forkWithPOR(k): # Check if the new quiver P is a fork with por k
                    self.couldBeFinite = False
                    continue

                if P.preForkWithPOR(k): # Check if the new quiver P is a pre-fork with por k
                    self.couldBeFiniteFP = False
                    self.couldBeFinite = False
                    continue
                    
                self.edges[Q][P] = k # If neither a fork or pre-fork with por k, then we can add it to our edges

                if P not in self.edges:
                    # Update potentional mutation invariants
                    self.couldBeMutationVortexFree = self.couldBeMutationVortexFree and not P.vortex()
                    self.couldBeMutationComplete = self.couldBeMutationComplete and not P.complete()
                    self.couldBeMutationCyclic = self.couldBeMutationCyclic and not P.acyclic()
                    self.mutationAcyclic = True if P.acyclic() else self.mutationAcyclic
                    
                    if self.couldBeMutationCyclic and self.mutationCyclicSubquiver is None:
                        if P.hasMutCyclicSubquiver():
                            self.foundCyclicSubquiver(P)
                    
                    # Update vertices, edges, and forefront
                    self.edges[P] = {Q : k}
                    self.vertices.append(P)
                    newForefront.append(P)
                    
                    # Update the representative lists
                    if P.numEdges <= self.leastEdges:
                        I = isomorphismClass(P, self.perms)[0]
                        if P.numEdges < self.leastEdges:
                            self.possibleReps = []
                            self.possibleIsoReps = []
                            self.leastEdges = P.numEdges

                        self.possibleReps.append(copy.deepcopy(P))
                        self.possibleIsoReps.append(I)
                else:
                    # This only occurs if the quiver we find is not a new one
                    self.edges[P][Q] = k

        self.forefront = newForefront # Update the forefront to be the new quivers we saw this round

        if len(self.forefront) == 0:
            self.finitePFP = True
            self.finite = self.couldBeFinite
            self.finiteFP = self.couldBeFiniteFP
            self.mutationAcyclic = not self.couldBeMutationCyclic
            self.hasVortex = not self.couldBeMutationVortexFree
            self.mutationComplete = self.couldBeMutationComplete

        self.size = len(self.vertices)

        return newForefront

    def updateInvariants(self):
        # Systematically updates all the invariants
        pass

    def foundCyclicSubquiver(self, Q):
        """Sets MutationCyclicSubquiver, not(MutationAcyclic), and sets the witness to Q."""
        self.mutationCyclicSubquiver = True
        self.mutationAcyclic = False
        self.mutationCyclicSubquiverWitness = Q


def isolatedQuiver(n):
    # returns the quiver with 0 arrows between n vertices
    
    m = [[0 for i in range(n)] for j in range(n)]

    return Quiver(m,False)

def generateLowWeightQuivers(n):
    # Generates the set of all possible quivers with low weights (|w| <= 2) of a given rank n

    seed = isolatedQuiver(n)
    
    result = [seed]

    lowWeights = [-2,-1,0,1,2]

    def oneStep(l, i, j):
        # takes in a list of quivers and produces a new list of quivers
        # where the (i,j) (index from 0) is replaced by every low weight
        # Does the same with (j,i) in a skew symmetric manner
        newList = []
        for Q in l:
            newM = copy.deepcopy(Q.matrix)
            for w in lowWeights:
                newM[i,j] = w
                newM[j,i] = -w
                newQ = Quiver(copy.deepcopy(newM))
                newList.append(newQ)

        return newList

    for i in range(n):
        for j in range(i+1,n):
            result = oneStep(result, i, j)

    return result

def permutations(n):
    # Gives a list of permutations on n
    #  proposal: itertools.permutations(range(n))
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    perm = permutations(n-1)
    newPerm = []
    for p in perm:
        first = [n-1] + p[:]
        last = p[:] + [n-1]
        middle = [p[:i] + [n-1] + p[i:] for i in range(1,n-1)]
        newP = [first]
        newP.extend(middle)
        newP.append(last)
        newPerm.extend(newP)

    return newPerm

def sign_perm(permutation):
    # Returns the sign of the given permutation.
    not_visited = set(permutation)
    cycle_sum = 0
    i = not_visited.pop()
    cycle = [i]
    while len(not_visited)>0:
        if permutation[i] in cycle:
            i = not_visited.pop()
            cycle_sum += len(cycle)-1
            cycle = [i]
        else:
            i = permutation[i]
            not_visited.remove(i)
            cycle.append(i)
    return (-1)**((cycle_sum + len(cycle)-1)% 2)


def isomorphicQuiver(quiver, p):
    # Returns the quiver given by the graph isomorphism perm
    
    n = quiver.n

    newMatrix = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(i+1,n):
            newMatrix[i][j] = quiver.matrix[p[i],p[j]]
            newMatrix[j][i] = quiver.matrix[p[j],p[i]]

    return Quiver(newMatrix, False)

def isomorphismClass(quiver, perms):
    # Gets the isomorphism class of a quiver
    # Returns a list of isomorphic quivers
    # Sorts list so that the greatest lexicographical quiver is at the end
    isoClass = [isomorphicQuiver(quiver, p) for p in perms]
    isoClass.sort()

    return isoClass

def reduceByIsomorphism(quivers):
    # Takes in a list of a quivers and returns a dictionary of all the distinct
    # isomorphism classes with the least lexi element as the representative
    # The value of the dict item is the number of times it appears in the list
    n = quivers[0].n
    perms = permutations(n)
    
    classes = [isomorphismClass(q, perms)[0] for q in quivers]

    return Counter(classes)

def removeDisconnected(quivers):
    # Removes all quivers which are not connected
    return [q for q in quivers if not q.connected()]

def removeAcyclic(quivers):
    # Removes all quivers that are acyclic
    return [q for q in quivers if not q.acyclic()]

def removeVortices(quivers):
    # Removes all quivers that are vortices
    return [q for q in quivers if not q.vortex()]

def boxQuiver(a,b):
    # Returns the box quiver with sides a and b
    M = [[0 for i in range(4)] for j in range(4)]

    for i in range(4):
        j = (i+1)% 4
        M[i][j] = a if i % 2 == 0 else b
        M[j][i] = - M[i][j]

    return Quiver(M)

def dreadedTorus():
    # Returns the box quiver with sides a and b
    M = [[] for j in range(4)]

    M[0] = [0,1,1,-1]
    M[1] = [-1,0,1,-1]
    M[2] = [-1,-1,0,2]
    M[3] = [1,1,-2,0]

    return Quiver(M)

def test():
    M = [[0,-2,-1,2],[2,0,0,-2],[1,0,0,0],[-2,2,0,0]]
    Q = Quiver(M)
    print(Q.determinant())
    classes = [mutationClass(P, permutations(4)) for P in isomorphismClass(Q, permutations(4))]

    for C in classes:
        for i in range(8):
            if not C.mutationCyclic:
                raise Exception(f"Happened here for {C} at mutation {i+1}")
            C.update()

    raise Exception("Testing finished")


def main():
    n = 4
    perms = permutations(n)
    box = isomorphismClass(boxQuiver(2,2),perms)[0]
    torus = isomorphismClass(dreadedTorus(),perms)[0]
    print("Generating low weight quivers:")
    quivers = generateLowWeightQuivers(n)
    numQuivers = len(quivers)
    print(f"{len(quivers)} low weight quivers generated. Reducing by isomorphism")
    reduced = reduceByIsomorphism(quivers)
    numIsoClasses = len(reduced)
    print(f"{len(reduced)} isomorphism classes remaining. Removing disconnected quivers:")
    reduced = [q for q in reduced if q not in removeDisconnected(reduced)]
    numConnected = len(reduced)
    print(f"{len(reduced)} connected quivers remaining. Removing Acyclic quivers:")
    reduced = removeAcyclic(reduced)
    numCyclic = len(reduced)
    print(f"{len(reduced)} quivers containing a cycle remaining. Removing quivers that are vortices:")
    reduced = removeVortices(reduced)
    numNonVortex = len(reduced)
    print(f"{len(reduced)} non-vortices remaining. Removing the box quiver and dreaded torus:")
    reduced.remove(box)
    reduced.remove(torus)
    numNonSpecial = len(reduced)
    print(f"{len(reduced)} quivers to test for mutation-acyclicity remaining.")
    m = 12
    print(f"Testing remaining by mutating up to {m} times.")
    mutClasses = [mutationClass(Q,perms) for Q in reduced]
    mutAcyclic = 0
    numVortex = 0
    numFinite = 0
    numFiniteFP = 0
    numFinitePFP = 0
    numMutCyclicSubquiver = 0
    mutAcyclicClasses = []
    mutCyclicSubquiverClasses = []
    mutationFiniteClasses = []
    mutationFiniteFPClasses = []
    mutationFinitePFPClasses = []
    mutationVortexClasses = []
    for M in mutClasses:
        for i in range(m):
            M.update()
            if not M.couldBeMutationCyclic:
                mutAcyclic += 1
                mutAcyclicClasses.append(M)
                break
            elif M.mutationCyclicSubquiver is not None and M.mutationCyclicSubquiver:
                numMutCyclicSubquiver += 1
                mutCyclicSubquiverClasses.append(M)
                break
                #print(M.rep)
                #print('---------')
            elif M.finite is not None and M.finite:
                numFinite += 1
                mutationFiniteClasses.append(M)
                break
            elif not M.couldBeMutationVortexFree:
                numVortex += 1
                mutationVortexClasses.append(M)
                #print(f"Found {numVortex} quivers with vortices in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finiteFP is not None and M.finiteFP:
                numFiniteFP += 1
                mutationFiniteFPClasses.append(M)
                #print(f"Found {numFiniteFP} quivers with Finite Forkless Part in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finitePFP is not None and M.finitePFP:
                numFinitePFP += 1
                mutationFinitePFPClasses.append(M)
                #print(f"Found {numFinitePFP} quivers with Finite Pre-Forkless Part in the mutation-class so far. Must be mutation-cyclic.")
                break

    numMutationCyclic = numNonSpecial - mutAcyclic
    specialClasses = mutAcyclicClasses + mutCyclicSubquiverClasses + mutationVortexClasses
    specialClasses = specialClasses + mutationFiniteClasses + mutationFiniteFPClasses + mutationFinitePFPClasses

    testing = [M for M in mutClasses if M not in specialClasses]

    print(f"Total: {numIsoClasses}")
    print(f"Connected: {numConnected} Disconnected: {numIsoClasses-numConnected}")
    print(f"Acyclic: {numConnected - numCyclic} Cyclic: {numCyclic}")
    print(f"Vortices: {numCyclic - numNonVortex} Non-vortices: {numNonVortex}")
    print(f"Special: 2 Non-special: {numNonSpecial}")
    print(f"Mutation-acyclic: {mutAcyclic} Mutation-cyclic: {numMutationCyclic}")
    print(f"Mutation Cyclic Subquiver: {numMutCyclicSubquiver}")
    print(f"Known Vortex Cyclic: {numVortex} Known Finite Cyclic: {numFinite}")
    print(f"Known FFP Cyclic: {numFiniteFP} Known FPFP Cyclic: {numFinitePFP}")
    print(f"Empirically Cyclic: {numMutationCyclic - numVortex - numFinite - numFiniteFP - numFinitePFP - numMutCyclicSubquiver}")

    print("Testing the Empirically Cyclic:")

    count = Counter(testing)
    print(f"There are {len(count)} many mutation classes without isomorphism")
    isoCount = Counter([M.isomorphicRepresentative() for M in testing])
    print(f"There are {len(isoCount)} many mutation classes with isomorphism")

    determinantCount = Counter([M.determinant for M in testing])
    numWithMutCyclicSubquiver = 0

    for M in isoCount:
        print(M)
        print("---------")
        print(M.determinant())
        print("---------")

    print("Number of distinct determinants: ", len(determinantCount))


if __name__ == "__main__":
    #test()
    main()