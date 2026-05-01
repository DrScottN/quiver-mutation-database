from Quiver import *
import sys
import time

LINEBREAK = 80*"#"

def reduceByIsomorphism(quivers):
    """
    Takes in a list of quivers and returns a list of all the distinct quivers up to isomorphism
    """
    # Takes in a list of a quivers and returns a dictionary of all the distinct
    # isomorphism classes with the least lexi element as the representative
    # The value of the dict item is the number of times it appears in the list
    n = quivers[0].n
    perms = permutations(n)
    
    classes = set(isomorphismClass(q, perms)[0] for q in quivers)

    return list(classes) # Counter(classes)

def betterReduceByIsomorphism(quivers):
    """
    Takes in a list of quivers and edits the list to remove all but the least isomorphic representatives
    """
    n = quivers[0].n
    perms = permutations(n)

    def onePass(perms, Q, quiverList):
        isoClass = isomorphismClass(Q, perms)
        leastQuiver = isoClass[0]

        return [P for P in quiverList if P not in isoClass], leastQuiver

    Q = quivers[0]
    newQuivers = []

    while len(quivers) > 0:
        quivers, P = onePass(perms, Q, quivers)
        newQuivers.append(P)
        Q = quivers[0]
        print(len(newQuivers))

    return newQuivers

def removeDisconnected(quivers):
    # Removes all quivers which are not connected
    return [q for q in quivers if q.connected()]

def removeAcyclic(quivers):
    # Removes all quivers that are acyclic
    return [q for q in quivers if not q.acyclic()]

def removeVortices(quivers):
    # Removes all quivers that have vortices
    return [q for q in quivers if q.vortex_free()]

def removeQuiversWithMutCyclicThreeCycle(quivers):
    # Removes all quivers that have a mutation-cyclic 3-cycle as a subquiver
    return [q for q in quivers if not q.hasMutCyclicThreeCycle()]

def specialQuivers(n, perms):
    quivers = []

    if n == 4:
        quivers.append(isomorphismClass(boxQuiver(2,2),perms)[0])
        quivers.append(isomorphismClass(dreadedTorus(),perms)[0])
    elif n == 5:
        pass
    else:
        print(f"No special quivers currently set for {n}")

    return quivers

def reduceQuivers(quivers, method):
    """
    Reduces the list of quivers by the method given
    Returns either a list or a Counter as well as the number of such quivers
    """
    reduced = method(quivers)

    return reduced, len(reduced)

def generateQuivers(n, weight_max):
    """
    Generates the whole class of quivers that we will be dealing with,
    Returned as a list as well as number of quivers generated
    """
    quivers = generateLowWeightQuivers(n, weight_max=weight_max)

    return quivers, len(quivers)

def isomorphismReductionFirst(n, perms, reduced):
    results = dict()
    print(f"{len(reduced)} low weight quivers generated. Reducing by isomorphism")
    reduced, numIsoClasses = reduceQuivers(reduced, reduceByIsomorphism)
    results["Isomorphism Classes"] = numIsoClasses

    print(f"{numIsoClasses} mutation classes up to isomorphism remaining. Removing disconnected quivers:")
    reduced, numConnected = reduceQuivers(reduced, removeDisconnected)
    results["Connected"] = numConnected
    results["Disconnected"] = numIsoClasses - numConnected

    print(f"{numConnected} connected quivers remaining. Removing Acyclic quivers:")
    reduced, numCyclic = reduceQuivers(reduced, removeAcyclic)
    results["Cyclic"] = numCyclic
    results["Acyclic"] = numConnected - numCyclic
   
    print(f"{numCyclic} quivers containing a cycle remaining. Removing quivers that contain a mutation-cyclic 3-cycle:")
    reduced, numNonMutCyclicThreeCycle = reduceQuivers(reduced, removeQuiversWithMutCyclicThreeCycle)
    results["Does Not Contain Mutation-cyclic Three Cycle"] = numNonMutCyclicThreeCycle
    results["Contains Mutation-cyclic Three Cycle"] = numCyclic - numNonMutCyclicThreeCycle

    print(f"{numNonMutCyclicThreeCycle} quivers that do not have a mutation-cyclic 3-cycle remaining. Removing quivers that have a vortex:")
    reduced, numNonVortex = reduceQuivers(reduced, removeVortices)
    results["Vortex-free"] = numNonVortex
    results["Contains a Vortex"] = numNonMutCyclicThreeCycle - numNonVortex
    
    print(f"{numNonVortex} non-vortices remaining. Removing the exceptional quivers:")
    
    for Q in specialQuivers(n, perms):
        if Q in reduced:
            reduced.remove(Q)
    
    numNonSpecial = len(reduced)
    results["Non-special Quivers"] = numNonSpecial
    results["Special Quivers"] = numNonVortex - numNonSpecial

    print(f"{numNonSpecial} quivers to test for mutation-acyclicity remaining.")
    
    print(LINEBREAK)

    return reduced, results

def isomorphismReductionLast(n, perms, reduced):
    results = dict()
    numQuivers = len(reduced)
    print(f"{numQuivers} low weight quivers generated. Removing disconnected quivers:")
    reduced, numConnected = reduceQuivers(reduced, removeDisconnected)
    results["Connected"] = numConnected
    results["Disconnected"] = numQuivers - numConnected

    print(f"{numConnected} connected quivers remaining. Removing Acyclic quivers:")
    reduced, numCyclic = reduceQuivers(reduced, removeAcyclic)
    results["Cyclic"] = numCyclic
    results["Acyclic"] = numConnected - numCyclic
   
    print(f"{numCyclic} quivers containing a cycle remaining. Removing quivers that contain a mutation-cyclic 3-cycle:")
    reduced, numNonMutCyclicThreeCycle = reduceQuivers(reduced, removeQuiversWithMutCyclicThreeCycle)
    results["Does Not Contain Mutation-cyclic Three Cycle"] = numNonMutCyclicThreeCycle
    results["Contains Mutation-cyclic Three Cycle"] = numCyclic - numNonMutCyclicThreeCycle

    print(f"{numNonMutCyclicThreeCycle} quivers that do not have a mutation-cyclic 3-cycle remaining. Removing quivers that have a vortex:")
    reduced, numNonVortex = reduceQuivers(reduced, removeVortices)
    results["Vortex-free"] = numNonVortex
    results["Contains a Vortex"] = numNonMutCyclicThreeCycle - numNonVortex
    
    print(f"{numNonVortex} non-vortices remaining. Removing the exceptional quivers:")
    
    for Q in specialQuivers(n, perms):
        if Q in reduced:
            reduced.remove(Q)
    
    numNonSpecial = len(reduced)
    results["Non-special Quivers"] = numNonSpecial
    results["Special Quivers"] = numNonVortex - numNonSpecial
    
    print(f"{numNonSpecial} quivers remaining. Reducing by isomorphism:")
    reduced, numIsoClasses = reduceQuivers(reduced, reduceByIsomorphism)
    results["Isomorphism Classes"] = numIsoClasses

    print(f"{numIsoClasses} quivers to test for mutation-acyclicity remaining.")
    
    print(LINEBREAK)

    return reduced, results



def findEmpiricallyCyclicMutClasses(reduced, perms, results, mutations):
    
    m = mutations
    print(f"Testing remaining by mutating up to {m} times.")
    print(LINEBREAK)
    
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
            if M.mutationAcyclic is True:
                mutAcyclic += 1
                mutAcyclicClasses.append(M)
                break
            elif M.mutationCyclicSubquiver is True:
                numMutCyclicSubquiver += 1
                mutCyclicSubquiverClasses.append(M)
                break
            elif M.finite is True:
                numFinite += 1
                mutationFiniteClasses.append(M)
                break
            elif M.hasVortex is True:
                numVortex += 1
                mutationVortexClasses.append(M)
                break
            elif M.finiteFP is True:
                numFiniteFP += 1
                mutationFiniteFPClasses.append(M)
                break
            elif M.finitePFP is True:
                numFinitePFP += 1
                mutationFinitePFPClasses.append(M)
                break

    numMutationCyclic = len(reduced) - mutAcyclic
    specialClasses = mutAcyclicClasses + mutCyclicSubquiverClasses + mutationVortexClasses
    specialClasses = specialClasses + mutationFiniteClasses + mutationFiniteFPClasses + mutationFinitePFPClasses

    results["Mutation-Acyclic"] = mutAcyclic
    results["Mutation-equivalent to a Quiver Containing a Mutation-cyclic Three Cycle"] = numMutCyclicSubquiver
    results["Mutation-Finite"] = numFinite
    results["Mutation-equivalent to a Quiver Containing a Vortex"] = numVortex
    results["Mutation-Finite Forkless Part"] = numFiniteFP
    results["Mutation-Finite Pre-forkless Part"] = numFinitePFP
    results["Emperically Mutation-Cyclic"] = len(reduced) - len(specialClasses)

    testing = [M for M in mutClasses if M not in specialClasses]

    return testing, results


def findDistinctMutClasses(testing):
    print(LINEBREAK)

    print("Testing the Empirically Cyclic:")
    print(LINEBREAK)

    count = Counter(testing)
    print(f"There are {len(count)} many mutation classes without isomorphism")
    isoCount = Counter([M.isomorphicRepresentative() for M in testing])
    print(f"There are {len(isoCount)} many mutation classes with isomorphism")

    determinantCount = Counter([M.determinant for M in testing])
    numWithMutCyclicSubquiver = 0

    print("Number of distinct determinants: ", len(determinantCount))

def alternateMain(n = 4, weight_max=2, mutations=5):
    perms = permutations(n)
    m = mutations
    
    print(LINEBREAK)
    print(f"Generating low weight quivers with {n} vertices and weights at most {weight_max}:")
    print(LINEBREAK)

    reduced, numQuivers = generateQuivers(n, weight_max)

    results = {"Total Quivers" : numQuivers}

    #reduced, reductionResults = isomorphismReductionFirst(n, perms, reduced)
    reduced, reductionResults = isomorphismReductionLast(n, perms, reduced)
    
    mutClasses, mutationResults = findEmpiricallyCyclicMutClasses(reduced, perms, reductionResults, m)

    results.update(reductionResults)
    results.update(mutationResults)
    for k, v in results.items():
        print(k, ":", str(v))

    findDistinctMutClasses(mutClasses)

def main(n=4, weight_max=2, mutations=5):
    perms = permutations(n)
    m = mutations
    
    print(LINEBREAK)
    print(f"Generating low weight quivers with {n} vertices and weights at most {weight_max}:")
    print(LINEBREAK)
    reduced, numQuivers = generateQuivers(n, weight_max)
    
    #print(f"{numQuivers} low weight quivers generated. Reducing by isomorphism")
    #reduced, numIsoClasses = reduceQuivers(reduced, reduceByIsomorphism)
    
    #print(f"{numIsoClasses} mutation classes up to isomorphism remaining. Removing disconnected quivers:")
    #reduced, numConnected = reduceQuivers(reduced, removeDisconnected)
    
    print(f"{numQuivers} low weight quivers generated. Removing disconnected quivers:")
    reduced, numConnected = reduceQuivers(reduced, removeDisconnected)

    print(f"{numConnected} connected quivers remaining. Removing Acyclic quivers:")
    reduced, numCyclic = reduceQuivers(reduced, removeAcyclic)
   
    print(f"{numCyclic} quivers containing a cycle remaining. Removing quivers that contain a mutation-cyclic 3-cycle:")
    reduced, numNonMutCyclicThreeCycle = reduceQuivers(reduced, removeQuiversWithMutCyclicThreeCycle)

    print(f"{numNonMutCyclicThreeCycle} quivers that do not have a mutation-cyclic 3-cycle remaining. Removing quivers that have a vortex:")
    reduced, numNonVortex = reduceQuivers(reduced, removeVortices)
    
    print(f"{numNonVortex} non-vortices remaining. Removing the exceptional quivers:")
    
    for Q in specialQuivers(n, perms):
        if Q in reduced:
            reduced.remove(Q)
    
    numNonSpecial = len(reduced)
    

    print(f"{numNonSpecial} non-exceptional quivers remaining. Reducing by isomorphism")
    reduced, numIsoClasses = reduceQuivers(reduced, reduceByIsomorphism)
    #print(f"{numNonSpecial} quivers to test for mutation-acyclicity remaining.")
    print(f"{numIsoClasses} isomorphism classes of quivers to test for mutation-acyclicity remaining.")
    print(LINEBREAK)
    

    print(f"Testing remaining by mutating up to {m} times.")
    print(LINEBREAK)
    
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
            if M.mutationAcyclic is True:
                mutAcyclic += 1
                mutAcyclicClasses.append(M)
                break
            elif M.mutationCyclicSubquiver is True:
                numMutCyclicSubquiver += 1
                mutCyclicSubquiverClasses.append(M)
                break
                #print(M.rep)
                #print('---------')
            elif M.finite is True:
                numFinite += 1
                mutationFiniteClasses.append(M)
                break
            elif M.hasVortex is True:
                numVortex += 1
                mutationVortexClasses.append(M)
                #print(f"Found {numVortex} quivers with vortices in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finiteFP is True:
                numFiniteFP += 1
                mutationFiniteFPClasses.append(M)
                #print(f"Found {numFiniteFP} quivers with Finite Forkless Part in the mutation-class so far. Must be mutation-cyclic.")
                break
            elif M.finitePFP is True:
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
    print(f"Mutation Cyclic 3-cycle: {numCyclic - numNonMutCyclicThreeCycle} Without Mutation Cyclic 3-cycle: {numNonMutCyclicThreeCycle}")
    print(f"Vortices: {numNonMutCyclicThreeCycle - numNonVortex} Non-vortices: {numNonVortex}")
    print(f"Special: {numNonVortex - numNonSpecial} Non-special: {numNonSpecial}")
    print(LINEBREAK)
    
    print(f"Now Listing Information Found by The {m} Mutations:")
    print(LINEBREAK)

    print(f"Mutation-acyclic: {mutAcyclic} Mutation-cyclic: {numMutationCyclic}")
    print(f"Mutation Cyclic Subquiver: {numMutCyclicSubquiver}")
    print(f"Known Vortex Cyclic: {numVortex} Known Finite Cyclic: {numFinite}")
    print(f"Known FFP Cyclic: {numFiniteFP} Known FPFP Cyclic: {numFinitePFP}")
    print(f"Empirically Cyclic: {numMutationCyclic - numVortex - numFinite - numFiniteFP - numFinitePFP - numMutCyclicSubquiver}")
    print(LINEBREAK)

    print("Testing the Empirically Cyclic:")
    print(LINEBREAK)

    count = Counter(testing)
    print(f"There are {len(count)} many mutation classes without isomorphism")
    isoCount = Counter([M.isomorphicRepresentative() for M in testing])
    print(f"There are {len(isoCount)} many mutation classes with isomorphism")

    determinantCount = Counter([M.determinant for M in testing])
    numWithMutCyclicSubquiver = 0

    #for M in isoCount:
    #    print(M)
    #    print("---------")
    #    print(M.determinant())
    #    print("---------")

    print("Number of distinct determinants: ", len(determinantCount))

if __name__ == "__main__":
    startTime = time.time()
    args = [int(j) for j  in sys.argv[1:]]
    if len(args) == 3:
        #main(*args)
        alternateMain(*args)
    else:
        main()
    endTime = time.time()
    print(LINEBREAK)
    print(f"Program took {endTime - startTime} seconds to run")
    print(LINEBREAK)
