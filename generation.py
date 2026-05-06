from quiver import *
import sys
import time
import csv

LINEBREAK = 80*"#"

def writeMutationClasses(muClassList, filename, indexLabel=False, initialQLabel=False):
    """ 
    Write the quivers from the given mutation class with their properties. 
    If indexLabel is True, include an entry for the mutation class itself (its index in the list muClassList)
    If initialQLabel is True, include the matrix of the first quiver is possibleReps with each entry. 
    """
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        header = headerForLabels()
        if indexLabel:
            header.append("class index")
        if initialQLabel:
            header.append("representative exchange matrix")
        writer.writerow(header)
        for i,C in enumerate(muClassList):
            for Qr in C.dataIterator():
                if indexLabel:
                    Qr.append(str(i))
                if initialQLabel:
                    Qr.append(str(C.possibleReps[0].matrix.flatten()))
                writer.writerow(Qr)

def readCSVDatabase(filename, header=True):
    """ 
    Iterator over given csv file.
    If header, treats the first row as dict labels. Automatically parsing any fields which agree with our labels
     and yield a dicitonary of the results.
    Otherwise, yield a list from each row, parsing only the first column as a quiver.
    """
    results = []
    with open(filename, 'r', newline='') as f:
        if header:
            reader = csv.DictReader(f)
            parseRow = parseMuClassHeaders
        else:
            reader = csv.reader(f)
            parseRow = lambda r : [parseSquareMatrixToQuiver(r[0])] + r[1:]
        for row in reader:
            yield parseRow(row)

def lowWeightQuiverIterator(n, maxWeight):
    """ Iterator for quivers on n vertices with small weights. """
    for w in itertools.product(range(-maxWeight, maxWeight+1), repeat=int(n*(n-1)/2)):
        yield Quiver([[(i<j)*w[n*i+j-((i+2)*(i+1)//2)] - (j<i)*w[n*j+i-((j+2)*(j+1)//2)] for j in range(n)] for i in range(n)])

def generateClasses(n, maxWeight, depth):
    """ 
    Generates all connected mutation classes with weights at most maxWeight, mutates them to depth, and gets the distinct unions. 
    """
    classes = []
    for Q in lowWeightQuiverIterator(n, maxWeight):
        if not Q.connected():
            continue
        CQ = mutationClass(Q, perms=permutations(n))
        for _d in range(depth):
            CQ.update()
        matched=False
        for i in range(len(classes)): #could keep these sorted and match faster. 
            if CQ == classes[i]:
                if matched:
                    classes[matched_index] = classes[matched_index].union(classes[i])
                    remove.append(i)
                    continue
                classes[i] = classes[i].union(CQ)
                matched=True
                matched_index=i
                remove = []
        if not matched:
            classes.append(CQ)
        else:
            for i in remove[::-1]:
                del classes[i]
    return classes

def generateClassesWithKnownAcyclicity(n, maxWeight, depth):
    """ Generates all connected mutation classes labeled with mutation-acyclicity (or not), with weights at most maxWeight, mutates them to depth. """
    for Q in lowWeightQuiverIterator(n, maxWeight):
        if not Q.connected():
            continue
        CQ = mutationClass(Q, perms=permutations(n))
        for _d in range(depth):
            CQ.update()
        if not (CQ.mutationAcyclic is Unknown):
            yield CQ


def reduceByIsomorphism(quivers):
    """
    Takes in a list of quivers and returns a list of all the distinct quivers up to isomorphism
    """
    n = quivers[0].n
    
    classes = set([isomorphismRep(q) for q in quivers])

    return list(classes)

def removeDisconnected(quivers):
    """
    Removes all quivers which are not connected
    """
    return [q for q in quivers if q.connected()]

def removeAcyclic(quivers):
    """
    Removes all quivers that are acyclic
    """
    return [q for q in quivers if not q.acyclic()]

def removeVortices(quivers):
    """
    Removes all quivers that have vortices
    """
    return [q for q in quivers if q.vortexFree()]

def removeQuiversWithMutCyclicThreeCycle(quivers):
    """
    Removes all quivers that have a mutation-cyclic 3-cycle as a subquiver
    """
    return [q for q in quivers if not q.hasMutCyclicThreeCycle()]

def specialQuivers(n, perms):
    """
    Produces all the 'special' quivers that we run across and needed to deal with. Currently not implemented for rank 5
    """
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

def generateQuivers(n, maxWeight):
    """
    Generates the whole class of quivers that we will be dealing with,
    Returned as a list as well as number of quivers generated
    """
    quivers = generateLowWeightQuivers(n, maxWeight=maxWeight)

    return quivers, len(quivers)

def isomorphismReductionFirst(n, perms, reduced):
    """
    Reduces the quivers given in reduced by first reducing by isomorphism
    It then steadily throws out quivers until the only thing we are left with
    Are quivers containing a cycle that we don't immediately know are mutation-cyclic
    
    Returns the reduced list and a results dictionary
    The results are keys of the properties with values of the number of quivers satisfying that property
    """

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
    """
    Reduces the quivers given in reduced by first steadily throwing out quivers until
    The only thing we are left with are quivers containing a cycle that we don't immediately
    Know are mutation-cyclic.
    It then reduces those quivers by isomorphism.
    
    Returns the reduced list and a results dictionary
    The results are keys of the properties with values of the number of quivers satisfying that property
    """
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
    """
    Goes through the reduced list of quivers and mutates all of them up to the mutations number
    It then returns the list of all mutation classes that we failed to show were mutation-acyclic or not
    These are then 'empirically' cyclic.

    Additionally, the results dictionary is added to and returned as well
    """
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
    """
    Takes in a list of mutation classes that are assumed to be empirically cyclic.
    It then forms a counter of these mutation classes, and finds the number of distinct
    Mutation classes up to mutation.
    It then does the same for the number of mutation classes up to mutation and isomorphism.
    Returns nothing
    """
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

def main(n = 4, maxWeight=2, mutations=5):
    """
    Function to find how many distinct mutation classes of quivers that we cannot 
    Prove mutation-acyclic or mutation-cyclic exist.
    Takes in a number of vertices, the max arrow weight, and the number of mutations to perform.
    Prints the results it finds along the way (Can be piped to file for saving)

    Returns nothing
    """

    perms = permutations(n)
    m = mutations
    
    print(LINEBREAK)
    print(f"Generating low weight quivers with {n} vertices and weights at most {maxWeight}:")
    print(LINEBREAK)

    reduced, numQuivers = generateQuivers(n, maxWeight)

    results = {"Total Quivers" : numQuivers}

    #reduced, reductionResults = isomorphismReductionFirst(n, perms, reduced)
    reduced, reductionResults = isomorphismReductionLast(n, perms, reduced) # Produces notable speedup
    
    mutClasses, mutationResults = findEmpiricallyCyclicMutClasses(reduced, perms, reductionResults, m)

    results.update(reductionResults)
    results.update(mutationResults)
    for k, v in results.items():
        print(k, ":", str(v))

    findDistinctMutClasses(mutClasses)


def generate2ClassDataset(acyclicQ, cyclicR, depth=3, labels=False, forklessDepth=3):
    """
    Iterator for a dataset of quivers mutation equivalent to acyclicQ and cyclicR out to depth.
    """
    classQ = mutationClass(acyclicQ, perms=permutations(acyclicQ.n))
    classR = mutationClass(cyclicR, perms=permutations(cyclicR.n))
    if labels:
        for d in range(forklessDepth): #compute invariants/data
            classQ.update()
            classR.update()
        labelsQ = classQ.labels()
        labelsR = classR.labels()
    #compute all quivers in range depth
    todo = [(acyclicQ,True,None,0), (cyclicR,False,None,0)] #quiver, label, last mu, depth
    while len(todo)>0:
        Q,l,p,d = todo.pop()
        for i in range(Q.n):
            if i==p:
                continue
            Qp = Q.mutate(i)
            if l: 
                yield (Qp, labelsQ) if labels else (Qp, l)
            else: 
                yield (Qp, labelsR) if labels else (Qp, l)
            if d < depth:
                todo.append((Qp,l,i,d+1))

def generateMulticlassDataset(classReps, depth, labels=False, forklessDepth=3):
    """ Iterator for mutation equivalence datasets with classes given by classReps. Assumes reps are mutation-distinct.
    If labels, then also computes the forkless part out to depth forklessDepth and uses the mutationclass labels. """
    if labels:
        classes = [mutationClass(Q, perms=permutations(Q.n)) for Q in classReps]
        for d in range(forklessDepth): #compute invariants/data
            for cl in classes:
                cl.update()
        labels = [cl.labels() for cl in classes]
    #compute all quivers in range depth
    todo = [(Q,i,None,0) for i,Q in enumerate(classReps)] #quiver, label, last mu, depth
    while len(todo)>0:
        Q,l,p,d = todo.pop()
        for i in range(Q.n):
            if i==p:
                continue
            Qp = Q.mutate(i)
            yield (Qp, labels[l]) if labels else (Qp, l)
            if d < depth:
                todo.append((Qp,l,i,d+1))

if __name__ == "__main__":
    startTime = time.time()
    args = [int(j) for j  in sys.argv[1:]]
    if len(args) == 3:
        main(*args)
    else:
        main()
    endTime = time.time()
    print(LINEBREAK)
    print(f"Program took {endTime - startTime} seconds to run")
    print(LINEBREAK)
