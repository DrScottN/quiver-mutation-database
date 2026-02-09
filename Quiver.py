import copy
from collections import Counter
import numpy as np 

class Quiver():
    def __init__(self, matrix, validate = True):
        # matrix is just a list of lists

        if validate:
            if len(matrix) != len(matrix[0]):
                raise Exception("Not a square matrix")
            else:
                for i, l in enumerate(matrix):
                    for j, a in enumerate(l):
                        if a != -matrix[j][i]:
                            raise Exception("Not a skew-symmetric matrix")

        self.matrix = matrix
        self.n = len(matrix)
        self.vertices = [v for v in range(self.n)]
        self.numEdges = sum(self.matrix[i][j] for i in self.vertices for j in self.vertices if self.matrix[i][j] > 0)

    def __hash__(self):
        t = tuple([self.matrix[i][j] for i in range(self.n) for j in range(i+1,self.n)])
        return hash(t)

    def __str__(self):
        return '\n'.join(str(r) for r in self.matrix)

    def __eq__(self, other):
        return self.matrix == other.matrix

    def __lt__(self, other):
        # Lexicographical order on matrix elements
        if self.__eq__(other):
            return False
        elif self.n != other.n:
            raise Exception("Attempted to compare two quivers of different size")

        for i in range(self.n):
            for j in range(self.n):
                if self.matrix[i][j] < other.matrix[i][j]:
                    return True
                elif self.matrix[i][j] > other.matrix[i][j]:
                    return False

        return False

    def __deepcopy__(self, memo):
        newMatrix = [[0 for i in range(self.n)] for j in range(self.n)]

        for i in range(self.n):
            for j in range(self.n):
                newMatrix[i][j] = self.matrix[i][j]

        return Quiver(newMatrix, validate=False)

    def determinant(self):
        if self.n % 2 == 1:
            return 0
        
        def pfaffian(matrix):
            n = len(matrix)

            if n == 0:
                return 1
            
            #print(n)
            return sum([pow(-1,j+1) * matrix[0][j] * pfaffian([[matrix[i][k] for i in range(1,n) if i != j] for k in range(1,n) if k != j]) for j in range(1,n)])

        return pfaffian(self.matrix)**2

    def updateWeight(self, w, i, j):
        # Updates the weight at index (i,j) and (j,i)
        # Returns the resulting modified quiver
        # Doesn't seem to modify the hash value though, which is concerning
        self.matrix[i][j] = w
        self.matrix[j][i] = -w
        return self
    
    def subquiverRemoveOneVertex(self, v):
        # Takes the subquiver by removing the vertex v
        return Quiver([[self.matrix[i][j] for i in range(self.n) if i != v] for j in range(self.n) if j != v])

    def connected(self):
        # Finds if the quiver is connected as a simple graph
        seen = [0] + [j for j in range(1,self.n) if self.matrix[0][j] != 0]
        numSeen = 1

        while numSeen < len(seen):
            extension = []
            for j in seen[numSeen:]:
                extension.extend([k for k in range(self.n) if self.matrix[j][k] != 0 and k not in extension])

            numSeen = len(seen)
            seen.extend([k for k in extension if k not in seen])

        return numSeen == self.n
    
    def sources(self):
        # Get the sources in a quiver

        return [i for i in range(self.n) if sum(self.matrix[i]) == sum([abs(self.matrix[i][j]) for j in range(self.n)])]
    
    def sinks(self):
        # Get the sources in a quiver

        return [i for i in range(self.n) if -sum(self.matrix[i]) == sum([abs(self.matrix[i][j]) for j in range(self.n)])]

    def acyclic(self):
        # Finds if hte quiver is acyclic or not
        if self.n < 3:
            return True

        sources = [i for i in range(self.n) if sum(self.matrix[i]) == sum([abs(self.matrix[i][j]) for j in range(self.n)])]

        if len(sources) == 0:
            return False
        elif len(sources) == self.n:
            return True
        
        vertices = [i for i in range(self.n) if i not in sources]

        newMatrix = [[self.matrix[i][j] for i in vertices] for j in vertices]

        return Quiver(newMatrix).acyclic()

    def complete(self):
        # Finds if the quiver is complete or not
        terms = [self.matrix[i][j] for i in range(self.n) for j in range(i+1,self.n)]

        s = 1

        for t in terms:
            s *= t

        return s != 0

    def abundant(self):
        # Finds if the quiver is abundant
        for i in range(self.n):
            for j in range(i+1, self.n):
                if abs(self.matrix[i][j]) < 2:
                    return False

        return True 
    
    def forkWithPOR(self, r):
        # Finds the quiver is a fork with point of return r
        if self.acyclic():
            return False
        elif not self.abundant(): 
            return False
        elif not self.subquiverRemoveOneVertex(r).acyclic():
            return False
        
        for i in range(self.n):
            if i == r:
                continue

            if self.matrix[i][r] > 0:
                for j in range(self.n):
                    if j == r or j == i:
                        continue

                    if self.matrix[r][j] > 0:
                        value = self.matrix[j][i] > 0
                        value = value and self.matrix[i][r] < self.matrix[j][i]
                        value = value and self.matrix[r][j] < self.matrix[j][i]

                        if not value:
                            return False
                        
            if self.matrix[i][r] < 0:
                for j in range(self.n):
                    if j == r or j == i:
                        continue

                    if self.matrix[r][j] < 0:
                        value = self.matrix[i][r] > self.matrix[j][i]
                        value = value and self.matrix[r][j] > self.matrix[j][i]
                        value = value and self.matrix[j][i] < 0

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

        if not Q.forkWithPOR(r) or not Q.forkWithPOR(r):
            return False

        for k in self.vertices:
            if self.matrix[i][k]*self.matrix[k][j] > 0:
                return False

        return True

    def preForkWithPOR(self, r):
        # Returns true if pre-fork with point of return r
        # else false (even if it is a fork with por r)
        
        for i in self.vertices:
            for j in self.vertices:
                if j == i or r in [i,j]:
                    continue

                if self.matrix[i][j] > 1:
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
            
            return (matrix[0][1]**2) + (matrix[1][2]**2) + (matrix[2][0]**2) - abs(matrix[0][1]*matrix[1][2]*matrix[2][0])

        return any(Q.threeCycle() and markovInvariant(Q.matrix) <= 4 and Q.abundant() for Q in [self.subquiverRemoveOneVertex(i) for i in self.vertices])

    def vortex(self):
        # Returns whether the quiver is a vortex or not
        if self.n != 4:
            return False
        
        vertices = self.sources() + self.sinks()

        if len(vertices) != 1:
            return False
        
        vertex = vertices[0]
        
        if 0 in self.matrix[vertex][:vertex] or 0 in self.matrix[vertex][vertex+1:]:
            return False
        
        Q = self.subquiverRemoveOneVertex(vertices[0])

        return Q.threeCycle()

    def oppositeQuiver(self):
        return Quiver([[-self.matrix[i][j] for i in range(self.n)] for j in range(self.n)])
    
    def mutate(self, k):
        # Gives the quiver formed by mutating at vertex k

        newMatrix = copy.deepcopy(self.matrix)

        for i in range(self.n):
            for j in range(self.n):
                if k in [i,j]:
                    newMatrix[i][j] = -newMatrix[i][j]
                elif self.matrix[i][k]*self.matrix[k][j] > 0:
                    add = self.matrix[i][k]*self.matrix[k][j] if self.matrix[i][k] > 0 else -self.matrix[i][k]*self.matrix[k][j]
                    newMatrix[i][j] += add

        return Quiver(newMatrix)

    def hasSourceSink(self):
        # Checks whether a source or sink exists

        return len(self.sources()) + len(self.sinks()) > 0     

    def chordless_cycles(self):
        # Iterator for all chordless cycles in the quiver. Formatted as list of vertex indices.
        
        def neighbors(v):
            return [j for j in range(self.n) if self.matrix[v][j] != 0]
        
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
                lefts += self.matrix[C[i]][C[(i+1)%len(C)]] < 0
            winding_dict[tuple(C)] = (wind - lefts, lefts)
        return winding_dict


    def cyclic_order(self):
        # Return the cyclic ordering of the quiver with all chordless cycles having winding number 1 (if oriented) or 0 (if acyclic). 
        #  Returns False if no such order exists. Only implemented for 4 vertex quivers.
        if self.n != 4:
            raise Exception("Not implemented for other than four vertices")
        if self.vortex():
            return False
        for sigma_p in permutations(3): #[[0,1,2,3], [0,2,1,3], [0,1,3,2], [0,2,3,1], [0,3,1,2], [0,3,2,1]]:
            sigma = sigma_p + [3]
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
        U = [[-(i<j)*self.matrix[sigma[i]][sigma[j]] for j in range(self.n)] for i in range(self.n)]
        for i in range(self.n):
            U[i][i]=1
        return U


    def markov(self):
        #computes the markov invariant. Return False if this quiver has no potentailly totally proper order.
        U = self.Umatrix()
        if U == False:
            return False
        #compute trace of U U^-1 (= U (I + N + N^2 + ...))
        raise Exception("Not implemented; no matrix inverse/multiplication yet.")

        

        


class mutationClass():
    def __init__(self, Q, perms, fast = False):
        # Takes in a quiver Q as the first member of our mutation class
        self.rep = copy.deepcopy(Q)
        self.vertices = [self.rep]
        self.forefront = [self.rep]
        self.edges = {self.rep : dict()} # self.edges[Q][P] gives the label of the mutation taking one to the other
        self.isoRep = isomorphismClass(self.rep, perms)[0]
        self.perms = perms
        self.possibleReps = [self.rep]
        self.possibleIsoReps = [self.isoRep]
        self.leastEdges = self.rep.numEdges
        self.complete = False
        self.acyclic = self.rep.acyclic()
        self.finiteFP = False
        self.couldBeFiniteFP = True
        self.finitePFP = False
        self.couldBeFinitePFP = True
        self.finite = False
        self.couldBeFinite = True
        self.hasVortex = self.rep.vortex()
        self.size = 1
        self.fast = fast
        self.determinant = self.rep.determinant()
        self.mutationComplete = self.rep.complete()
        self.mutationCyclic = self.rep.hasMutCyclicSubquiver()

    def __str__(self):
        return str(self.representative())

    def __hash__(self):
        return hash(self.representative())

    def __eq__(self, other):
        return self.representative() == other.representative()

    def isomorphic(self, other):
        return self.isomorphicRepresentative() == other.isomorphicRepresentative()

    def representative(self):
        if len(self.possibleReps) > 1:
            self.possibleReps.sort()
        
        return self.possibleReps[0]

    def isomorphicRepresentative(self):
        if len(self.possibleIsoReps) > 1:
            self.possibleIsoReps.sort()

        return self.possibleIsoReps[0]

    def update(self):
        # updates the mutation class by mutating every quiver on the edge in every possible direction
        # Naturally throws away forks and pre-forks if they have the same point of return
        if self.complete:
            return None

        newForefront = []

        for Q in self.forefront:
            #print(Q)
            #print('--------')
            for k in Q.vertices:
                P = Q.mutate(k)
                #print("Mutating at k = ", k)
                #print(P)
                #print('-------')
                if not P.forkWithPOR(k):
                    self.edges[Q][P] = k # this will cause an error, but hopefully not an important one
                    if P.preForkWithPOR(k):
                        self.couldBeFiniteFP = False
                        self.couldBeFinite = False
                        continue

                    if P not in self.edges:
                        self.acyclic = self.acyclic or P.acyclic()
                        self.hasVortex = self.hasVortex or P.vortex() 
                        self.edges[P] = {Q : k}
                        if not self.fast:
                            self.vertices.append(P)
                            self.mutationComplete = self.mutationComplete or P.complete()
                            self.mutationCyclic = self.mutationCyclic or P.hasMutCyclicSubquiver()

                        
                        newForefront.append(P)
                        
                        if P.numEdges <= self.leastEdges:
                            I = isomorphismClass(P, self.perms)[0]
                            if P.numEdges < self.leastEdges:
                                self.possibleReps = []
                                self.possibleIsoReps = []
                                self.leastEdges = P.numEdges

                            self.possibleReps.append(copy.deepcopy(P))
                            self.possibleIsoReps.append(I)
                    else:
                        self.edges[P][Q] = k
                else:
                    self.couldBeFinite = False


        self.forefront = newForefront

        if len(self.forefront) == 0:
            self.complete = True
            self.finitePFP = True
            self.finite = self.couldBeFinite
            self.finiteFP = self.couldBeFiniteFP
            self.mutationCyclic = not self.acyclic

        self.size = len(self.vertices)

        return newForefront

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
                newM[i][j] = w
                newM[j][i] = -w
                newQ = Quiver(copy.deepcopy(newM))
                newList.append(newQ)

        return newList

    for i in range(n):
        for j in range(i+1,n):
            result = oneStep(result, i, j)

    return result

def permutations(n):
    # Gives a list of permutations on n
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

def isomorphicQuiver(quiver, p):
    # Returns the quiver given by the graph isomorphism perm
    
    n = quiver.n

    newMatrix = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(i+1,n):
            newMatrix[i][j] = quiver.matrix[p[i]][p[j]]
            newMatrix[j][i] = quiver.matrix[p[j]][p[i]]

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

if __name__ == "__main__":
    #test()
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
    for M in mutClasses:
        for i in range(m):
            M.update()
            if M.acyclic:
                mutAcyclic += 1
                #print(f"Found {mutAcyclic} mutation-acyclic quivers so far.")
                break
            elif M.mutationCyclic:
                numMutCyclicSubquiver += 1
                break
                #print(M.rep)
                #print('---------')
            elif M.finite:
                numFinite += 1
                break
            elif M.hasVortex:
                numVortex += 1
                #print(f"Found {numVortex} quivers with vortices in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finiteFP:
                numFiniteFP += 1
                #print(f"Found {numFiniteFP} quivers with Finite Forkless Part in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finitePFP:
                numFinitePFP += 1
                #print(f"Found {numFinitePFP} quivers with Finite Pre-Forkless Part in the mutation-class so far. Must be mutation-cyclic.")
                break

    numMutationCyclic = numNonSpecial - mutAcyclic

    testing = [M for M in mutClasses if not(M.acyclic or M.hasVortex or M.finitePFP or M.mutationCyclic or M.finite or M.finitePFP)]

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
