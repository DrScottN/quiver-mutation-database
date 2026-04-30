from Quiver import *

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

def main(n=4, weight_max=2, mutations=5):
    #n = 4
    perms = permutations(n)
    box = isomorphismClass(boxQuiver(2,2),perms)[0]
    torus = isomorphismClass(dreadedTorus(),perms)[0]
    print("Generating low weight quivers:")
    quivers = generateLowWeightQuivers(n, weight_max=weight_max)
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
    if box in reduced: reduced.remove(box)
    if torus in reduced: reduced.remove(torus)
    numNonSpecial = len(reduced)
    print(f"{len(reduced)} quivers to test for mutation-acyclicity remaining.")
    m = mutations
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
    main()
