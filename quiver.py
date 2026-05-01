import copy
from collections import Counter
import numpy as np 
import itertools
import math

import polynomial
from unknown import unknown

#custom entry for 'none, but we might find it later!'
Unknown = unknown()


class Quiver():
    def __init__(self, matrix, validate = True):
        if validate and not isinstance(matrix, np.ndarray):
            if len(matrix) != len(matrix[0]):
                raise Exception("Not a square matrix")
            else:
                for i, l in enumerate(matrix):
                    for j, a in enumerate(l):
                        if a != -matrix[j][i]:
                            raise Exception("Not a skew-symmetric matrix")

        self.matrix = np.array(matrix, dtype=object) # type=object lets us put arbitrary python in there.
        self.n = len(matrix)
        self.vertices = [v for v in range(self.n)]
        self.numEdges = sum(self.matrix[i,j] for i in self.vertices for j in self.vertices if self.matrix[i,j] > 0)

    def __hash__(self):
        t = tuple(self.matrix[i,j] for i in range(self.n) for j in range(i+1,self.n))
        return hash(t)

    def __str__(self):
        return '\n'.join(str(r) for r in self.matrix)

    def __eq__(self, other):
        return np.array_equal(self.matrix, other.matrix) #.all()

    def __lt__(self, other):
        """ 
        Lexicographical order on matrix elements
        assumes skew-symmetric matrix
        """
        if self.n != other.n:
            raise Exception("Attempted to compare two quivers of different size")

        for i in range(self.n):
            for j in range(i+1,self.n):
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
            
            return sum([pow(-1,j+1) * matrixAsList[0][j] * pfaffian([[matrixAsList[i][k] for i in range(1,n) if i != j] for k in range(1,n) if k != j]) for j in range(1,n)])

        return pfaffian(self.matrix.tolist())**2
    
    def subquiver(self, Vs):
        """ 
        Gives the subquiver of self with vertex set Vs
        Vs is an enumeration of elements from 1 to n, vertex labels are not preserved.
        currently no rule against Vs=[0,0,0], which would give an isolated quiver;
        result is only a subquiver if Vs is a set of distinct elements. 
        """
        return Quiver([[self.matrix[i,j] for i in Vs] for j in Vs]) 

    def subquiverRemoveOneVertex(self, v):
        """ Takes the subquiver by removing the vertex v """
        return self.subquiver([i for i in range(self.n) if i != v]) # Quiver([[self.matrix[i,j] for i in range(self.n) if i != v] for j in range(self.n) if j != v])


    def connected(self):
        """ Finds if the quiver is connected as a simple graph """
        seen = [0]
        #unseen = set([])
        for j in range(1, self.n):
            if self.matrix[0,j] != 0:
                seen.append(j)
            #else:
                #unseen.add(j)
        #seen = [0] + [j for j in range(1,self.n) if self.matrix[0,j] != 0]
        numSeen = 1
        while numSeen < len(seen):
            go_to = len(seen)
            for j in range(numSeen,go_to):
                for k in range(1,self.n):
                    if self.matrix[seen[j],k] != 0 and k not in seen:
                        seen.append(k)

            numSeen = go_to
            if len(seen) == self.n:
                return True

        return False #len(seen) == self.n
    
    def sources(self):
        """ Get the sources in a quiver """
        return [i for i in range(self.n) if all([self.matrix[i,j]>=0 for j in range(self.n)])]
    
    def sinks(self):
        """ Get the sinks in a quiver """

        return [i for i in range(self.n) if all([self.matrix[i,j]<=0 for j in range(self.n)])]

    def acyclic(self):
        """ Finds if the quiver is acyclic or not """
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
        """ Finds if the quiver is complete or not """
        terms = [self.matrix[i,j] for i in range(self.n) for j in range(i+1,self.n)]

        s = 1

        for t in terms:
            s *= t

        return s != 0

    def abundant(self):
        """ Finds if the quiver is abundant """
        for i in range(self.n):
            for j in range(i+1, self.n):
                if abs(self.matrix[i,j]) < 2:
                    return False

        return True 
    
    def forkWithPOR(self, r):
        """ returns true if a fork with point of return r """
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
        """ returns true if pre-fork with given vertices """
        if self.acyclic():
            return False
        
        Q = self.subquiverRemoveOneVertex(i)
        P = self.subquiverRemoveOneVertex(j)

        if not Q.forkWithPOR(r -(i<r)) or not P.forkWithPOR(r - (j<r)):
            return False

        for k in self.vertices:
            if self.matrix[i,k]*self.matrix[k,j] > 0:
                return False

        return True

    def preForkWithPOR(self, r):
        """
        Returns true if pre-fork with point of return r
        else false (even if it is a fork with por r)
        """
        for i,j in itertools.combinations([x for x in self.vertices if x != r], 2):
            if self.matrix[i,j] > 1:
                continue
            if self.preForkWithVertices(r, i, j):
                return True
        return False

    def threeCycle(self):
        """ Returns whether the quiver is a 3-cycle or not """
        if self.n != 3:
            return False
        
        return not self.acyclic()

    def hasMutCyclicThreeCycle(self):
        """ Returns whether the quiver has a 3-cycle as a subquiver """
        def markovInvariant(matrix):
            if len(matrix) != 3:
                raise Exception("Not implemented yet")
            
            return (matrix[0,1]**2) + (matrix[1,2]**2) + (matrix[2,0]**2) - abs(matrix[0,1]*matrix[1,2]*matrix[2,0])

        return any(Q.threeCycle() and markovInvariant(Q.matrix) <= 4 and Q.abundant() for Q in [self.subquiver(s) for s in itertools.combinations(self.vertices, 3)])

    def vortex(self):
        """ 
        Returns whether the quiver is a vortex or not.
         Vortices are defined in "Long Mutation Cycles" by FN
        """
        if self.n != 4:
            return False
        
        apexes = self.sources() + self.sinks()

        if len(apexes) != 1:
            return False
        apex = apexes[0]
        
        if 0 in self.matrix[apex,:apex] or 0 in self.matrix[apex,apex+1:]:
            return False
        Q = self.subquiverRemoveOneVertex(apex)
        return Q.threeCycle()

    def vortexFree(self):
        """ determines if the quiver has no vortex subquiver """
        for V in itertools.combinations(self.vertices, 4):
            if self.subquiver(V).vortex():
                return False
        return True

    def oppositeQuiver(self):
        return Quiver(-1*self.matrix)
    
    def mutate(self, k):
        """ Gives the quiver formed by mutating at vertex k """

        newMatrix = copy.deepcopy(self.matrix)

        for i in range(self.n):
            for j in range(self.n):
                if k in [i,j]:
                    newMatrix[i,j] = -newMatrix[i,j]
                elif self.matrix[i,k]*self.matrix[k,j] > 0:
                    add = self.matrix[i,k]*self.matrix[k,j] if self.matrix[i,k] > 0 else -self.matrix[i,k]*self.matrix[k,j]
                    newMatrix[i,j] += add
        return Quiver(newMatrix)

    def matrixMutate(self, k):
        """ Gives the quiver formed by mutating at vertex k via matrix mutation """
        # empirically slower than the mutation method above.
        C = np.identity(self.n, dtype='object')
        C[k,:]= np.vectorize(lambda x : max(x,0))(self.matrix[k,:])
        C[k,k]=-1
        return Quiver(np.transpose(C) @ self.matrix @ C)


    def hasSourceSink(self):
        """ Checks whether a source or sink exists """
        return len(self.sources()) + len(self.sinks()) > 0     

    def chordlessCycles(self):
        """ Iterator for all chordless cycles in the quiver. Formatted as list of vertex indices. """
        
        def neighbors(v):
            return [j for j in range(self.n) if self.matrix[v,j] != 0]
        
        def chordlessCyclesAt(v):
            """ Iterator for chordless cycles with minimal element v """
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
            for C in chordlessCyclesAt(i):
                yield tuple(C)

    def windingData(self, sigma):
        """ Using linear order sigma (a permutation of the indices), return winding numbers and edges against direction of traversal of the chordless cycles. """
        winding_dict = dict()
        for C in self.chordlessCycles():
            wind = 0 #times orbited
            lefts = 0 #arrows against flow
            for i in range(len(C)):
                wind += (sigma[C[(i+1)%len(C)]] < sigma[C[i]])
                lefts += self.matrix[C[i],C[(i+1)%len(C)]] < 0
            winding_dict[tuple(C)] = (wind - lefts, lefts)
        return winding_dict


    def cyclicOrder(self):
        """ 
        Return the cyclic ordering of the quiver with all chordless cycles having winding number 1 (if oriented) or 0 (if acyclic). 
         Returns False if no such order exists. 
        """
        if self.n==1:
            return [0]
        if self.n==2:
            return [0,1]
        for sigma_p in permutations(self.n-1): #[[0,1,2,3], [0,2,1,3], [0,1,3,2], [0,2,3,1], [0,3,1,2], [0,3,2,1]]:
            sigma = sigma_p + (self.n-1,)
            W = self.windingData(sigma)
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

    def proper(self):
        """ returns a proper cyclic ordering (but not necessarily totally proper) if it exists. """
        for p in permutations(self.n):
            winding = self.windingData(p)
            next_p = False
            for C, x in winding.items():
                ell = x[1]
                r = len(C)-ell
                if (ell > 1 and x[0] == 1-ell) or (r > 1 and x[0] == r-1):
                    next_p = True
                    break
            if next_p:
                continue
            return p
        return False


    def Umatrix(self):
        """ 
        Compute a unipotent companion. 
         Return False if this quiver has no potentailly totally proper order.
        """
        
        sigma = self.cyclicOrder()
        if sigma==False:
            return False
        U = np.array([[-(i<j)*self.matrix[sigma.index(i),sigma.index(j)] for j in range(self.n)] for i in range(self.n)], dtype=object)
        for i in range(self.n):
            U[i,i]=1
        return U


    def markov(self):
        """computes the markov invariant. 
         Return False if this quiver has no potentially totally proper order.
        """
        U = self.Umatrix()
        if U is False:
            return False
        #compute trace of U U^-T (= U (I + N + N^2 + ...))
        N = np.transpose(np.identity(self.n, dtype='object') - U)
        Cosquare = U.copy()
        for e in range(1,self.n):
            Cosquare += U @ np.linalg.matrix_power(N, e)
        trace = 0
        for i in range(self.n):
            trace += Cosquare[i,i]
        return self.n - trace
        
    def alexanderPolynomial(self):
        """ 
        Computes the alexander polynomial (expressed as a polynomial object). 
        Returns False if the quiver is not potentially totally proper.
        """
        U = self.Umatrix()
        if U is False:
            return False
        Ux = np.array([[U[i,j]*polynomial.polynomial([0,1]) for i in range(self.n)] for j in range(self.n)], dtype=object) - np.array([[U[j,i]*polynomial.polynomial([1]) for i in range(self.n)] for j in range(self.n)], dtype=object)
        det = polynomial.polynomial([0])
        for p in permutations(self.n):
            a = polynomial.polynomial([1])
            for i in range(self.n):
                a *= Ux[i, p[i]]
            det += a * signOfPerm(p)
        return det

    def gcdVector(self):
        """ Computes the gcd vector of self, as a tuple with ith entry the gcd of row i. """
        return tuple(math.gcd(*[self.matrix[i,j] for j in range(self.n)]) for i in range(self.n))

    def gcdTwoPath(self):
        """Computes the gcd vector of all pairs of (not-necessarily-oriented) paths."""
        d = 0
        for i in range(self.n):
            for j in range(i+1, self.n):
                for k in range(k, self.n):
                    t = math.gcd([abs(self.matrix[i,j]*self.matrix[j,k]), abs(self.matrix[j,i]*self.matrix[i,k]), abs(self.matrix[j,k]*self.matrix[k,i])])
                    d = math.gcd(d, t)
        return d
        
    def casalsDeterminant(self):
        """ return the mod-4 determinant of a quasi-cartan companion. This is a mutation invariant. """
        #  see arxiv: 2311.03601
        #  uses slow det formula.
        A = np.vectorize(abs)(self.matrix)
        for i in range(self.n):
            A[i,i] = 2
        det = 0
        for p in permutations(self.n):
            a = 1
            for i in range(self.n):
                a *= (A[i, p[i]])%4
            det = (det + a * signOfPerm(p))%4
        return det

    def sevenCongruence(self):
        """ return the 'corank' of A mod 4 up to congruence, via the dimension of a particular space """
        #  see https://www.sciencedirect.com/science/article/abs/pii/S0022404925000593 thm 1.7 (note it is written differently on arxiv)
        #  very slow, we enumerate vectors.
        V0 = []
        V000 = []
        A = np.vectorize(abs)(self.matrix)
        for i in range(self.n):
            A[i,i] = 2
        for v in itertools.product([0,1], repeat=self.n):
            v = np.array(v)
            wp = v@A
            if (wp%2==0).all():
                V0.append(v)
        for v in V0:
            for delta in itertools.product([0,2], repeat=self.n):
                vp = (v + np.array(delta))%4
                wp = vp@A
                if (wp%4==0).all():
                    V000.append(v)
                    break #try next v, we found a good delta.
        #return the dimension of V000 over F_2
        return len(V000).bit_length() - 1




def eigenvalues(M):
    """ 
    Computes the eigenvalues of the given matrix, cast to int64.
      Note that the rank is a mutation invariant but not these values. 
    """
    return np.linalg.eig(np.int64(M)).eigenvalues

def matrixRank(M):
    """ Computes the rank of a given matrix, cast to int64 """
    return np.linalg.matrix_rank(np.int64(M))

def descendFork(Q, maxSteps=-1):
    """If Q is a fork, mutate at the point of return until either we reach a non-fork or determine there is no forkless part.
    Returns the resulting quiver.
    If maxSteps is set, will apply at most maxSteps mutations."""
    prevMutation = -1 #the last mutation performed, to prevent stutter/loops.
    while maxSteps != 0:
        foundDescent = False
        maxSteps -= 1
        for r in range(Q.n):
            if r==prevMutation:
                continue
            if Q.forkWithPOR(r):
                Q = Q.mutate(r)
                foundDescent = True
                prevMutation = r
                break
        if not foundDescent:
            break
    return Q
    

def sinkSet(quiver):
    """ Calculates all sink/source mutation equivalent quivers """
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

def generateAcyclicsFromAlexander(alex):
    """ Takes an Alexander polynomial alex and generates all acyclic quivers with that polynomial.
     Iterator. Assumes alex is uni-variate. Not Fast. """
    if alex is False: return
    degree = max([j[0] for j in alex.coeffDict.keys() if alex.coeffDict[j]!=0])
    markov = degree + alex.coeffDict[(degree-1,)]
    if markov < 0 or markov is False:
        return

    for w in itertools.product(list(range(math.isqrt(markov)+1)), repeat=(degree*(degree-1))//2):
        # iterates over all combinations of weights < sqrt(markov), which is much too large
        #  faster is to build subquivers that are smaller than markov. Even faster potentially is to use more coefficients and a solver.
        
        X = np.array([[0 for i in range(degree)] for j in range(degree)])
        #print(degree, np.tril_indices(X.shape[0], k = -1), w)
        X[np.tril_indices(X.shape[0], k = -1)] = w
        Q = Quiver(X - np.transpose(X))
        if Q.alexanderPolynomial() != alex:
            continue
        #to avoid repeats; this is very slow.
        if any([any([R < Q for R in isomorphismClass(Rp, permutations(degree)) if (np.tril(R.matrix) >= 0).all()]) for Rp in sinkSet(Q)]):
            continue

        yield Q

def generateAcyclicsBelowMarkov(markov, n, connected=True):
    """Iterator for all acyclic quivers on n vertices with markov invariant <= markov.
    checks connected if connected=True."""
    for w in itertools.product(list(range(math.isqrt(markov)+1)), repeat=(n*(n-1))//2):
        X = np.zeros((n,n), dtype='object')
        X[np.tril_indices(n, k = -1)] = w
        Q = Quiver(X - np.transpose(X))
        if connected and not Q.connected():
            continue
        if Q.markov() <= markov:
            #to avoid repeats; this is very slow.
            if any([any([R < Q for R in isomorphismClass(Rp, permutations(degree)) if (np.tril(R.matrix) >= 0).all()]) for Rp in sinkSet(Q)]):
                continue
            yield Q


def surfaceQuiver(quiver):
    """ Return if the quiver is a surface quiver by identifying a block decomposition. """
    #  Block decomposition described in section 4 of "Skew-symmetric cluster algebras of finite mutation type" arxiv:0811.1703
    #  and in FST 13.1

    # easy check
    if (2<quiver.matrix).any():
        return False
    # constants
    n = quiver.n
    neighbors_out = [[j for j in range(n) if quiver.matrix[i,j] >0] for i in range(n)]
    neighbors_in = [[j for j in range(n) if quiver.matrix[i,j] <0] for i in range(n)]
    neighbor_count = [len(neighbors_in[i] + neighbors_out[i]) for i in range(n)]

    # more easy checks
    if max(neighbor_count) > 8:
        return False
    if max([len(neighbors_in[i]) for i in range(n)]) > 4:
        return False
    if max([len(neighbors_out[i]) for i in range(n)]) > 4:
        return False
    # constructors for the blocks. 
    #  B_ gives the matrix of the block, 
    #  v_ gives the covering count of each vertex in Q by that matrix 
    #   (v[i]==0 means i is not in any block, v[i]==2 means i is dead).

    #white vert
    def v0(i):
        return np.array([x==i for x in range(n)], dtype=np.int8)
    def B0(i):
        B = np.zeros((n,n), dtype='object')
        return B
    #single arrow, out out i->j
    def vI(i,l):
        return np.array([x in [i,l] for x in range(n)], dtype=np.int8)
    def BI(i,l):
        B = np.zeros((n,n), dtype='object')
        B[i,l] += 1
        return B - np.transpose(B)
    #oriented cycle, out out out i->j->k->i
    def vII(i,j,k):
        return np.array([x in [i,j,k] for x in range(n)], dtype=np.int8)
    def BII(i,j,k):
        B = np.zeros((n,n), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[k,i] += 1
        return B - np.transpose(B)
    #sink, out dead dead j->i<-k
    def vIII(i,j,k):
        return np.array([(x==i) + 2*(x in [j,k]) for x in range(n)], dtype=np.int8)
    def BIIIa(i,j,k):
        B = np.zeros((n,n), dtype='object')
        B[j,i] += 1
        B[k,i] += 1
        return B - np.transpose(B)
    #source, out dead dead j<-i->k
    def BIIIb(i,j,k):
        B = np.zeros((n,n), dtype='object')
        B[i,j] += 1
        B[i,k] += 1
        return B - np.transpose(B)
    #diamond, out dead dead out i->j->l->i->k->l
    def vIV(i,j,k,l):
        return np.array([(x in [i,l]) + 2*(x in [j,k]) for x in range(n)], dtype=np.int8)
    def BIV(i,j,k,l):
        B = np.zeros((n,n), dtype='object')
        B[i,j] += 1
        B[i,k] += 1
        B[j,l] += 1
        B[k,l] += 1
        B[l,i] += 1
        return B - np.transpose(B)
    #surrounded, out dead dead dead dead, i->j->k<-g->h<-j, k->i->g, i<-h
    def vV(i,j,k,g,h):
        return np.array([(x==i) + 2*(x in [j,k,g,h]) for x in range(n)], dtype=np.int8)
    def BV(i,j,k,g,h):
        B = np.zeros((n,n), dtype='object')
        B[i,j] += 1
        B[j,k] += 1
        B[j,h] += 1
        B[k,i] += 1
        B[h,i] += 1
        B[g,k] += 1
        B[g,h] += 1
        B[i,g] += 1
        return B - np.transpose(B)

    #init vars; these are shared among the recursive functions below via 'nonlocal'.
    B = np.zeros((n,n), dtype='object')
    v = np.zeros(n,dtype=np.int8)
    blocks = []
    def check_status():
        nonlocal B, v
        # check if B,v is consistent so far
        matched_verts = [i for i in range(n) if v[i]==2]
        # do finished vertices match?
        for i,j in itertools.combinations(matched_verts, 2):
            if B[i,j] != quiver.matrix[i,j]:
                return False
        return True
    
    def block_decomp():
        nonlocal B,v,blocks
        # recursive function for finding a block decomposition
        #  B is matrix so far, v counts how matched vertices are (0 is unused, 1 is outlet, 2 is dead)
        outlet_verts = []
        matched_verts = []
        unused_verts = []
        for i in range(n):
            if v[i]==1:
                outlet_verts.append(i)
            elif v[i]==2:
                matched_verts.append(i)
            else:
                unused_verts.append(i)
        if (B == quiver.matrix).all():
            return True

        def try_block(v_fun, B_fun, block):
            nonlocal B, v, blocks
            vp = v_fun(*block[1])
            if ((v+vp)>2).any():
                return False
            Bp = B_fun(*block[1])
            B = B+Bp
            v = v+vp
            if check_status():
                blocks.append(block)
            
                if block_decomp():
                    return True
                blocks.pop()
            B = B-Bp
            v = v-vp
            return False

        if len(outlet_verts)==0:
            u = unused_verts[0]
            if neighbor_count[u]==0: #isolated vertex
                v[u]=2
                return block_decomp()
            #attach one dead and recurse, or drop down
            #dead => BIII (1 nbr), BIV (2 nbr), BV (3 nbr)
            #BIII
            match neighbor_count[u]:
                case 1:
                    # u -> i <- j or u <- i -> j
                    i = (neighbors_in[u]+neighbors_out[u])[0]
                    if v[i]!=2:
                        if quiver.matrix[u,i] >0:
                            for j in [x for x in neighbors_in[i] if neighbor_count[x]==1 and v[x]==0 and x!=u]:
                                if try_block(vIII, BIIIa, ('IIIa',(i,u,j))): return True
                        else:
                            for j in [x for x in neighbors_out[i] if neighbor_count[x]==1 and v[x]==0 and x!=u]:
                                if try_block(vIII, BIIIb, ('IIIb',(i,u,j))): return True
                case 2:
                    if len(neighbors_in[u]) == len(neighbors_out[u]):
                        # i -> u -> l -> i -> k -> l
                        i,l = neighbors_in[u][0], neighbors_out[u][0]
                        if v[i]!=2 and v[l]!=2 and i in neighbors_out[l]:
                            for k in [x for x in neighbors_out[i] if x in neighbors_in[l] and v[x]==0]:
                                if try_block(vIV, BIV, ('IV', (i,u,k,l))): return True
                case 3: 
                    if len(neighbors_out[u])==1: 
                        #i->j->u<-g->h<-j
                        i = neighbors_out[u]
                        j,g = neighbors_in[u]
                        if len(neighbors_out[j])==2:
                            h = neighbors_out[j][0] if neighbors_out[j][1]==u else neighbors_out[j][1]
                            if try_block(vV, BV, ('V', (i,j,u,g,h))): return True
                        
                    elif len(neighbors_in[u])==1: 
                        #i->j->k<-u->h<-j
                        i = neighbors_in[u]
                        k,h = neighbors_out[u]
                        if len(neighbors_in[k])==2:
                            j = neighbors_in[k][0] if neighbors_in[k][1]==u else neighbors_in[k][1]
                            if try_block(vV, BV, ('V', (i,j,k,u,h))): return True
        else:
            u = outlet_verts[0]
        #now have u with at most one block attached.
        # outlet => B1 (i or j) B2 BIIIab BIV BV
        # this version tries with no assumptions; could instead check nbrs+current count
        for i in outlet_verts+unused_verts: #2+ outlets
            if i==u: continue
            if neighbor_count[u] < 6 and neighbor_count[i] < 6: #degree at most 5
                if i not in neighbors_out[u] and try_block(vI, BI, ('I',(i,u))): return True
                if i not in neighbors_in[u] and try_block(vI, BI, ('I',(u,i))): return True
            if neighbor_count[i]<7 and neighbor_count[u]<7: #degree at most 6
                for j in outlet_verts+unused_verts:
                    if j in [i,u]: continue
                    if j not in neighbors_out[u] and try_block(vII, BII, ('II',(u,i,j))): return True
                    if j not in neighbors_in[u] and try_block(vII, BII, ('II',(u,j,i))): return True
            if i not in neighbors_in[u]: #u->i
                for j,k in itertools.combinations([x for x in neighbors_out[i] if x in neighbors_in[u]], 2):
                    if try_block(vIV, BIV, ('IV',(i,j,k,u))): return True
            if i not in neighbors_out[u]: #i->u
                for j,k in itertools.combinations([x for x in neighbors_in[i] if x in neighbors_out[u]], 2):
                    if try_block(vIV, BIV, ('IV',(u,j,k,i))): return True
        for i,j in itertools.combinations([x for x in neighbors_out[u] if v[x]==0], 2): #1 outlet, outset
            match (neighbor_count[i], neighbor_count[j]):
                case (1,1):
                    if try_block(vIII, BIIIb, ('IIIb',(u,i,j))): return True
                case (3,3):
                    if len(neighbors_out[i])==2:
                        g,h = neighbors_out[i]
                        if try_block(vV, BV, ('V',(u,i,g,j,h))): return True
        for i,j in itertools.combinations([x for x in neighbors_in[u] if v[x]==0], 2):
            if try_block(vIII, BIIIa, ('IIIa',(u,i,j))): return True
        if v[u]==1: #try the skip block
            if try_block(v0, B0, ('0', (u,))): return True
        return

    if block_decomp():
        return [True, blocks]
    else:
        return False

class mutationClass():
    """ 
    Object to hold mutation-equivalent quivers
    Each quiver is represented as a vertex of the mutation class
    Each mutation is represented as an edge of the mutation class
    Edges are given as a dict of dicts. 
    For example, for a quiver Q one mutation away from P at vertex k:
      {Q : {P : k}, P : {Q : k}} 
    """
    # Could easily use this to build a graph visualization of mutation classes

    # The class will be split up into two main ideas: essentially a set for mutations and a tool for exploring mutation-classes

    def __init__(self, Q = None, perms = None, vertices = None, edges = None, maxWeight=2**16):
        """ 
        Takes in a quiver Q as the first member of our mutation class. perms is required.
         set maxWeight=False to ignore weight checks 
         (which stop mutations in a direction if a quiver with more than maxWeight arrows between any pair of vertices is encountered).
        """
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
        self.possibleIsoReps = [isomorphismRep(self.initialQ)] #[isomorphismClass(self.initialQ, perms)[0]]

        self.maxWeight = maxWeight
        self.hitMaxWeight = False
        for Q in self.vertices:
            if self.maxWeight and np.max(np.vectorize(abs)(Q.matrix)) > self.maxWeight:
                self.hitMaxWeight = True
                #do we wish to do anything here?


        # Now, check for the mutation-invariants we can immediately find
        self.mutationConnected = self.initialQ.connected()
        self.determinant = self.initialQ.determinant()
        self.leastEdges = self.initialQ.numEdges
        self.size = 1
        self.gcdVector = self.initialQ.gcdVector()
        self.BRank = matrixRank(self.initialQ.matrix)
        self.alexanderPolynomial = self.initialQ.alexanderPolynomial()
        self.casalsDeterminant = self.initialQ.casalsDeterminant()
        self.sevensCongruence = self.initialQ.sevenCongruence()

        if self.initialQ.acyclic():
            self.mutationAcyclic = True
            self.mutationCyclicSubquiver = False
            self.mutationCyclicSubquiverWitness = None
            self.totallyProper = True
        else:
            if self.initialQ.hasMutCyclicThreeCycle():
                self.foundCyclicSubquiver(self.initialQ)
            else:
                self.mutationCyclicSubquiver = Unknown
                self.mutationCyclicSubquiverWitness = None
                self.mutationAcyclic = Unknown
            self.totallyProper = Unknown if self.alexanderPolynomial is not False else False

        self.mutationComplete = False if self.initialQ.complete() is not True else Unknown
        self.isSurfaceQuiver = Unknown

        self.hasVortex = False if self.mutationAcyclic is True else Unknown
        self.finite = Unknown
        self.finiteFP = Unknown
        self.finitePFP = Unknown
        # Setup possible mutation-invariants
        for Q in self.vertices:
            if self.maxWeight and np.max(np.vectorize(abs)(Q.matrix)) > 2:
                self.isSurfaceQuiver = False
                self.finite = False
            
            if self.hasVortex is Unknown and not Q.vortexFree():
                self.hasVortex = True
                self.mutationAcyclic = False
        if self.isSurfaceQuiver is Unknown:
            self.isSurfaceQuiver = bool(surfaceQuiver(self.initialQ))
            self.finite = True if self.isSurfaceQuiver else self.finite
            self.finiteFP = True if self.isSurfaceQuiver else self.finiteFP
            self.finitePFP = True if self.isSurfaceQuiver else self.finitePFP

        self.hasMutationCycle = Unknown
        self.mutationAbundant = False if not self.initialQ.abundant() else Unknown

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
        """ 
        Allows us to initialize a mutationClass without starting with a single quiver
        We assume that vertices lists all vertices
        Edges may not have the edges for both sides
        Edges does not include vertices not present in vertices
        """

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

        if not self._consistentProperties(other):
            return None
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
        
        #we may know more about this than usual, update
        newMutClass._setProperties(self)
        newMutClass._setProperties(other)

        return newMutClass

    def intersection(self, *others):
        """ 
        Gets the intersection of two mutation classes
        Returns a mutation class containing the common quivers
        Else returns none if no quivers are common 
        """
        other = others[0]
        if len(others) > 1:
            other = other.intersection(*others[1:])

        if other is None:
            return None
        elif isinstance(other, Quiver):
            other = mutationClass(other, self.perms)

        if not self._consistentProperties(other):
            return None

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

        #we know more about this than usual, update
        newMutClass._setProperties(self)
        newMutClass._setProperties(other)
        
        return newMutClass

    def _setProperties(self, other):
        """ 
        Helper function that updates properties of this mutation class to match those in the other class.
         Assumes that the two classes intersect. Otherwise this can cause unreachable state. 
        """
        self.finite = self.finite if self.finite is not Unknown else other.finite
        self.finitePFP = self.finitePFP if self.finitePFP is not Unknown else other.finitePFP
        self.finiteFP = self.finiteFP if self.finiteFP is not Unknown else other.finiteFP
        self.hasVortex = self.hasVortex if self.hasVortex is not Unknown else other.hasVortex
        self.hasMutationCycle = self.hasMutationCycle if self.hasMutationCycle is not Unknown else other.hasMutationCycle
        self.mutationAbundant = self.mutationAbundant if self.mutationAbundant is not Unknown else other.mutationAbundant
        self.mutationAcyclic = self.mutationAcyclic if self.mutationAcyclic is not Unknown else other.mutationAcyclic
        self.totallyProper = self.totallyProper if self.totallyProper is not Unknown else other.totallyProper
        self.mutationComplete = self.mutationComplete if self.mutationComplete is not Unknown else other.mutationComplete
        self.mutationCyclicSubquiver = self.mutationCyclicSubquiver if self.mutationCyclicSubquiver is not Unknown else other.mutationCyclicSubquiver
        # handled by init: self.mutationConnected, determinant, gcdVector, BRank, alexanderPolynomial, casalsDeterminant, sevensCongruence, 

    def _consistentProperties(self, other):
        """check if self and other might intersect/union. helper function to skip pointless intersections/unions."""
        if self.mutationConnected != other.mutationConnected: return False
        if self.determinant != other.determinant: return False
        if self.gcdVector != other.gcdVector: return False
        if self.BRank != other.BRank: return False
        if self.casalsDeterminant != other.casalsDeterminant: return False
        if self.sevensCongruence != other.sevensCongruence: return False
        #doesn't help: if not consistent_ternary(self.finite, other.finite): return False
        return True

    def agrees(self, other):
        """Determine if self and other could be mutation equivalent, given current data."""
        if not self._consistentProperties(other):
            return False
        # check conditional invariants

        def consistentTernary(X,Y):
           """Determines if X and Y could agree"""
           if (X is Unknown) or (Y is Unknown):
               return True
           if X == Y:
               return True
           return False

        if not consistentTernary(self.totallyProper, other.totallyProper):
            False
        if (self.totallyProper | other.totallyProper) is True:
            if self.alexanderPolynomial != other.alexanderPolynomial:
                return False

        if not consistentTernary(self.mutationAcyclic, other.mutationAcyclic):
            False
        if not consistentTernary(self.mutationCyclicSubquiver, other.mutationCyclicSubquiver):
            False
        if not consistentTernary(self.mutationComplete, other.mutationComplete):
            False
        if not consistentTernary(self.isSurfaceQuiver, other.isSurfaceQuiver):
            False
        if not consistentTernary(self.hasVortex, other.hasVortex):
            False
        if not consistentTernary(self.finite, other.finite):
            False
        if not consistentTernary(self.finiteFP, other.finiteFP):
            False
        if not consistentTernary(self.finitePFP, other.finitePFP):
            False
        return True
        
        
    def labels(self):
        """
        generate tab seperated string of invariant data associated to this class
        """
        def encoding(T):
            if T is Unknown:
                return 'U'
            if T is True:
                return 'T'
            if T is False:
                return 'F'
            else:
                return str(T)
        data = []
        data.append(encoding(self.mutationAcyclic))
        data.append(encoding(self.mutationAbundant))
        data.append(encoding(self.mutationComplete))
        data.append(encoding(self.mutationConnected))
        data.append(encoding(self.mutationCyclicSubquiver))
        data.append(encoding(self.mutationFiniteClasses))
        data.append(encoding(self.mutationFiniteFPClasses))
        data.append(encoding(self.mutationFinitePFPClasses))
        data.append(encoding(self.isSurfaceQuiver))
        data.append(encoding(self.hasVortex))

        data.append(encoding(self.totallyProper))
        data.append(encoding(self.alexanderPolynomial))

        data.append(encoding(self.initialQ.n))
        data.append(encoding(self.determinant))
        data.append(encoding(self.gcdVector))
        data.append(encoding(self.BRank))
        data.append(encoding(self.casalsDeterminant))
        data.append(encoding(self.sevensCongruence))

        data.append(encoding(self.hasMutationCycle))
        data.append(encoding(self.leastEdges))
        return "\t".join(data)
        

    def emptyIntersection(self,other):
        return self.intersection(other) is None

    def representative(self):
        # Gets a quiver with lowest edge weights to represent the class
        if len(self.possibleReps) > 1:
            self.possibleReps.sort()
        
        return self.possibleReps[0]

    def isomorphicRepresentative(self):
        """ Gets an isomorphic quiver with lowest edge weights to represent the class """
        if len(self.possibleIsoReps) > 1:
            self.possibleIsoReps.sort()

        return self.possibleIsoReps[0]

    def update(self):
        """ 
        updates the mutation class by mutating every quiver on the edge in every possible direction
        Automatically throws away forks and pre-forks if they have the same point of return as where we mutated
        """        
        #This function has a lot of room to benefit from parallelization
        #As the quivers in the forefront are effectively independent from each other
        newForefront = []

        for Q in self.forefront:
            for k in Q.vertices:
                P = Q.mutate(k) # Mutate the Quiver
                if P.forkWithPOR(k): # Check if the new quiver P is a fork with por k
                    self.finite = False
                    continue

                if P.preForkWithPOR(k): # Check if the new quiver P is a pre-fork with por k
                    self.finiteFP = False
                    self.finite = False
                    continue
                    
                
                if self.maxWeight and np.max(np.vectorize(abs)(P.matrix)) > self.maxWeight:
                    self.hitMaxWeight = True
                    continue

                self.edges[Q][P] = k # If neither a fork or pre-fork with por k, then we can add it to our edges

                if P not in self.edges:
                    # Update potentional mutation invariants
                    self.hasVortex = True if not P.vortexFree() else self.hasVortex
                    self.mutationComplete = False if not P.complete() else self.mutationComplete
                    self.mutationAbundant = False if not P.abundant() else self.mutationAbundant
                    if self.mutationAcyclic is Unknown:
                        if P.acyclic():
                            self.mutationAcyclic = True
                            self.hasVortex=False
                            self.mutationAbundant = P.abundant()
                            self.totallyProper = True
                    
                    if (self.mutationAcyclic is Unknown or self.mutationAcyclic is False) and (self.mutationCyclicSubquiver is Unknown):
                        if P.hasMutCyclicThreeCycle():
                            self.foundCyclicSubquiver(P)
                    
                    # Update vertices, edges, and forefront
                    self.edges[P] = {Q : k}
                    self.vertices.append(P)
                    newForefront.append(P)
                    
                    # Update the representative lists
                    if P.numEdges <= self.leastEdges:
                        I = isomorphismRep(P) #isomorphismClass(P, self.perms)[0]
                        if P.numEdges < self.leastEdges:
                            self.possibleReps = []
                            self.possibleIsoReps = []
                            self.leastEdges = P.numEdges

                        self.possibleReps.append(copy.deepcopy(P))
                        self.possibleIsoReps.append(I)
                else:
                    # This only occurs if the quiver we find is not a new one
                    # Hence a mutation cycle has appeared
                    self.edges[P][Q] = k
                    self.hasMutationCycle = True


        self.forefront = newForefront # Update the forefront to be the new quivers we saw this round

        if len(self.forefront) == 0 and not self.hitMaxWeight: #if we pruned, don't be overconfident.
            self.finitePFP = True
            self.finite = True if self.finite is Unknown else self.finite
            self.finiteFP = True if self.finiteFP is Unknown else self.finiteFP
            self.mutationAcyclic = False if self.mutationAcyclic is Unknown else self.mutationAcyclic
            self.hasVortex = True if self.hasVortex is Unknown else self.hasVortex
            self.mutationComplete = False if self.mutationComplete is Unknown else self.mutationComplete
            self.hasMutationCycle = False if self.hasMutationCycle is Unknown else self.hasMutationCycle
            self.mutationAbundant = True if self.mutationAbundant is Unknown else self.mutationAbundant
            if not self.hasVortex and self.mutationComplete:
                self.totallyProper=True

        self.size = len(self.vertices)

        return newForefront

    def foundCyclicSubquiver(self, Q):
        """Sets MutationCyclicSubquiver, not(MutationAcyclic), and sets the witness to Q."""
        self.mutationCyclicSubquiver = True
        self.mutationAcyclic = False
        self.mutationCyclicSubquiverWitness = Q


def isolatedQuiver(n):
    """ returns the quiver with 0 arrows between n vertices """
    m = [[0 for i in range(n)] for j in range(n)]

    return Quiver(m,False)

def generateLowWeightQuivers(n, maxWeight=2):
    """ Generates the set of all possible quivers with low weights (|w| <= maxWeight) of a given rank n """
    seed = isolatedQuiver(n)
    
    result = [seed]

    lowWeights = list(range(-maxWeight, maxWeight+1)) #[-2,-1,0,1,2]

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
    """ Gives a list of permutations on n """
    return list(itertools.permutations(range(n)))

def signOfPerm(permutation):
    """ Returns the sign of the given permutation. """
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
    """ Returns the quiver given by the graph isomorphism perm """
    
    n = quiver.n

    newMatrix = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(i+1,n):
            newMatrix[i][j] = quiver.matrix[p[i],p[j]]
            newMatrix[j][i] = quiver.matrix[p[j],p[i]]

    return Quiver(newMatrix, False)

def isomorphismClass(quiver, perms):
    """ 
    Gets the isomorphism class of a quiver
    Returns a list of isomorphic quivers
    Sorts list so that the greatest lexicographical quiver is at the end
    """
    isoClass = [isomorphicQuiver(quiver, p) for p in perms]
    isoClass.sort()

    return isoClass

def isomorphismRep(Q):
    """Returns the lex min rep of Q in its isomorphism class"""
    #assumes perms is all the permutations.
    maxWeight = 0
    opts1=[]
    opts2=[]
    for i,j in itertools.combinations(range(Q.n), 2):
        v = abs(Q.matrix[i,j])
        if maxWeight < v:
            maxWeight = v
            opts1=[]
            opts2=[]
        if maxWeight == v:
            if Q.matrix[i,j] > 0: 
                opts1.append(j)
                opts2.append(i)
            else: 
                opts1.append(i)
                opts2.append(j)
    isoCandidates = []
    for i in range(len(opts1)):
        p1,p2 = opts1[i], opts2[i]
        L = [j for j in range(Q.n) if j not in [p1,p2]]
        for pTail in itertools.permutations(L):
            p = (p1, p2) + pTail
            isoCandidates.append(isomorphicQuiver(Q, p))
    return min(isoCandidates)
    

        
    



def boxQuiver(a,b):
    """ Returns the box quiver with sides a and b """
    M = [[0 for i in range(4)] for j in range(4)]

    for i in range(4):
        j = (i+1)% 4
        M[i][j] = a if i % 2 == 0 else b
        M[j][i] = - M[i][j]

    return Quiver(M)

def dreadedTorus():
    """ Returns the box quiver with sides a and b """
    M = [[] for j in range(4)]

    M[0] = [0,1,1,-1]
    M[1] = [-1,0,1,-1]
    M[2] = [-1,-1,0,2]
    M[3] = [1,1,-2,0]

    return Quiver(M)


def main():
    print("You should be using generation.py")

if __name__ == "__main__":
    #test()
    main()
