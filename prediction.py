from quiver import *
import random

# Contains algorithms and helper functions mirroring how we would solve mutation equivalence or mutation acyclicity in practice.
#  includes variations with more/less compute, certainty, and data. 

def invariantsDictionary(Q, includeSurface=True):
    """ return a dict of mutation invariants of Q.  """
    invs = dict()
    invs["rank"] = Q.n
    invs["connected"] = Q.connected()
    invs["B rank"] = matrixRank(Q.matrix)
    invs["Casal's Det"] = Q.casalsDeterminant()
    invs["gcd"] = Q.gcdVector() 
    if includeSurface: invs["surf quiver"] = bool(surfaceQuiver(Q))
    #note, block decomp is not an invariant but being surface type is
    return invs

def mutationEquivalentLocal(Q,R, upToIsomorphism=False):
    """ 
    return what we can determine about the mutation-equivalence of Q and R
    Does not mutate either quiver. 
     if upToIsomorphism is True, compares Q against all quivers isomorphic to R 
    """
    if upToIsomorphism:
        all_known = True
        for p in permutations(R.n):
            result = mutationEquivalentLocal(Q, isomorphicQuiver(R, p))
            if result is True:
                return result
            if result is Unknown:
                all_known = False
        if all_known:
            return False
        return Unknown
    
    if invariantsDictionary(Q) != invariantsDictionary(R):
        return False
    if Q.acyclic() or R.acyclic():
        if Q.alexanderPolynomial() != R.alexanderPolynomial():
            return False
    return Unknown


def mutationAcyclicLocal(Q, matchAlexander=False):
    """ 
    return what we can determine about the mutation-acyclicity of Q
    Does not mutate Q, move to forkless part first for best results.
     if matchAlexander=True, then generates all acyclic quivers with the same alexander poly
     and checks if they could be mutation equivalent (based on local data).
      This can be extremely slow for large quivers.
    suggested convention: round Unknown to matchAlexander 
    """
    if Q.acyclic():
        return True
    if Q.hasMutCyclicThreeCycle():
        return False
    if Q.markov() is False: #includes Q contains a vortex
        return False
    if Q.markov() < 0 or (Q.connected() and Q.markov() < Q.n-1):
        return False
    if Q.determinant()**2 > (2*Q.markov())**Q.n:
        return False
    if matchAlexander:
        for R in generateAcyclicsFromAlexander(Q.alexanderPolynomial()):
            if mutationEquivalentLocal(Q,R) is not False:
                return Unknown #Heuristically, likely to be True
        return False
    else:
        return Unknown

def mutationEquivalent(Q, R, distance, upToIsomorphism=False, returnClasses=False, Qclass=None, Rclass=None):
    """
    Determines if Q is isomorphic to R using ~distance mutations, along with the invariants.
    returns True if they are equivalent, False if they are not, and Unknown if we cannot determine either way.
    
    upToIsomorphism=True considers the problem up to isomorphism
    returnClasses=True returns a tuple (Result, Mutation_class_of_Q, Mutation_class_of_R)
    Q_class, R_class may be provided for efficiency; distance is ignored if a class is given.
    attempts to replace Q,R with non-forks in distance mutations if class is not provided. 
     This may result in more than distance mutations being applied.
    """
    if upToIsomorphism:
        permsQ = permutations(Q.n)
        permsR = permutations(R.n)
    else:
        permsQ = [list(range(Q.n))]
        permsR = [list(range(R.n))]

    update = False
    if Qclass is None:
        Q = descendFork(Q, maxSteps=distance)
        Qclass = mutationClass(Q, perms=permsQ)
        update = True
    if Rclass is None:
        R = descendFork(R, maxSteps=distance)
        Rclass = mutationClass(R, perms=permsR)
        update = True

    local_check = mutationEquivalentLocal(Q, R, upToIsomorphism=upToIsomorphism)
    if local_check is not Unknown:
        return (local_check, Qclass, Rclass) if returnClasses else local_check
    
    #we may now assume same rank, etc.

    if update:
        distance = (distance+1)//2
        for d in range(distance):
            R_Class.update()
            Q_class.update()
    
    #check if overlap found
    if union := Rclass.union(Qclass):
        return (True, union, union) if returnClasses else True

    #check if invariants are consistent
    if Rclass.agrees(Qclass) is False: #not implemented
        return (False, Qclass, Rclass) if returnClasses else False

    return (Unknown, Qclass, Rclass) if returnClasses else Unknown




def runBaselineAcyclic(dataset, depthSearch=1, train=None, useSurface=True, useExhaustive=False, seed=2842026):
    """Runs several variants of our baseline against the given dataset, and prints the accuracies achieved.
        dataset: iterable of tuples, (quiver, correct_label), quiver is a quiver object and correct_label is a bool
        train (optional): iterable of (quiver, correct_label) pairs 
        (train may be used to guess other labels, if None then some baselines will sample examples from dataset)
        useSurface: if True, computes the block decomposition for invariants.
        useExhaustive: if True, enumerates all acyclic quivers with relevant size.
        seed: used for random in selecting subset of data, irrelevant if train is not None.
        If the local check ever fails to classify a large number of quivers, this can be very slow."""

    dataset_copy = list(dataset) #could use itertool.tee(dataset, k) and run in parallel for memory savings
    D = len(dataset_copy)
    dataset_not_local = []
    acyclic_count = 0
    cyclic_count = 0
    random.seed(seed)

    correct=0
    U_T_correct=0
    U_F_correct=0
    # baseline 1: blind attempts. 
    print("Local tests")
    #  3 ways to round: unknown is incorrect; unknown as True, unknown takes False.
    for Q,l in dataset_copy:
        result = mutationAcyclicLocal(Q, matchAlexander=useExhaustive)
        if result is not Unknown:
            assert result == l, f"Dataset and baseline disagree; {Q.matrix} has label {l} but we find {result}."
            correct += 1
        else:
            dataset_not_local.append((Q,l))
            if l is True:
                U_T_correct+=1
            if l is False:
                U_F_correct+=1
        if l: 
            acyclic_count += 1
        else: 
            cyclic_count += 1

    print(f" Unknown->incorrect: {correct} ({100*correct/D}%)")
    print(f"*Unknown->True: {correct+U_T_correct} ({100*(correct+U_T_correct)/D}%)")
    print(f" Unknown->False: {correct+U_F_correct} ({100*(correct+U_F_correct)/D}%)")
    
    # baseline 2: with training data.
    print()
    print("Comparing against a train set.")
    if train is None:
        # select a random selection of quivers of each label
        #  we choose to take samples of ~1/20 of the whole dataset; a bit weak if dataset contains few mutations from many classes
        #  provide a comparable train set to your model otherwise.
        acyclic_egs = random.choices([d[0] for d in dataset_copy if d[1]], k=round(acyclic_count/20))
        cyclic_egs = random.choices([d[0] for d in dataset_copy if not d[1]], k=round(cyclic_count/20))
    else:
        acyclic_egs = [d[0] for d in train if d[1]]
        cyclic_egs = [d[0] for d in train if not d[1]]
    # build useful invariant vectors for each class
    acyclic_invariants = []
    for Q in acyclic_egs:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in acyclic_invariants:
            acyclic_invariants.append(d)

    cyclic_invariants = []
    for Q in cyclic_egs:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in cyclic_invariants:
            cyclic_invariants.append(d)

    print(f"Found {len(cyclic_invariants)} distinct cyclic invariants and {len(acyclic_invariants)} distinct acyclic invariants.")
    overlap = len([d for d in acyclic_invariants if d in cyclic_invariants])
    if overlap > 0:
        print(f"They intersect for {overlap} values.")
    else:
        print("The sets of invariants are mutually distinct.")
    print()

    #attempt to label the unlabeled with the invariants
    correct_inv_with_cyclic_Alexander = 0
    correct_inv_without_cyclic_Alexander = 0
    undetermined_with_cyclic_Alexander = [0,0] #pair, [#correct label is False, #correct label is True]
    undetermined_without_cyclic_Alexander = [0,0]
    for Q,l in dataset_not_local:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        #match if invariants are in acyclic, in cyclic (without alexander considered), and in cyclic with alexander considered
        match (d in acyclic_invariants, any([all([d[k] == inv[k] for k in d.keys()]) for inv in cyclic_invariants]), d in cyclic_invariants):
            case (True, True, True):
                #ambiguous still
                undetermined_with_cyclic_Alexander[l] += 1 #hack; True=1 as an int while False=0, so this increments the count according to label.
                undetermined_without_cyclic_Alexander[l] += 1
            case (True, False, False):
                if l==True:
                    correct_inv_without_cyclic_Alexander += 1
                    correct_inv_with_cyclic_Alexander += 1
            case (True, True, False):
                undetermined_without_cyclic_Alexander[l] += 1
                if l==True:
                    correct_inv_with_cyclic_Alexander += 1
            case (False, True, True):
                if l==False:
                    correct_inv_without_cyclic_Alexander += 1
                    correct_inv_with_cyclic_Alexander += 1
            case (False, True, False):
                undetermined_with_cyclic_Alexander[l] += 1
                if l==False:
                    correct_inv_without_cyclic_Alexander += 1
            case (False, False, False):
                undetermined_with_cyclic_Alexander[l] += 1
                undetermined_without_cyclic_Alexander[l] += 1
            case (_b, False, True):
                assert False, "Implementation error in heuristic."

    print("Ignoring Alexander Polynomial for cyclic examples:")
    c = correct + correct_inv_without_cyclic_Alexander
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f"*Unknown->True: {c+undetermined_without_cyclic_Alexander[True]} ({100*(c+undetermined_without_cyclic_Alexander[True])/D}%)")
    print(f" Unknown->False: {c+undetermined_without_cyclic_Alexander[False]} ({100*(c+undetermined_without_cyclic_Alexander[False])/D}%)")
    
    print("Using Alexander Polynomials for the cyclic examples:")
    c = correct + correct_inv_with_cyclic_Alexander
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f" Unknown->True: {c+undetermined_with_cyclic_Alexander[True]} ({100*(c+undetermined_with_cyclic_Alexander[True])/D}%)")
    print(f" Unknown->False: {c+undetermined_with_cyclic_Alexander[False]} ({100*(c+undetermined_with_cyclic_Alexander[False])/D}%)")

    # baseline 3: Apply mutations
    print()
    print("Applying mutations to our examples and checking for equivalence.")
    acyclic_eg_classes = [mutationClass(Q, perms=permutations(Q.n)) for Q in acyclic_egs]
    cyclic_eg_classes = [mutationClass(Q, perms=permutations(Q.n)) for Q in cyclic_egs]
    for QC in acyclic_eg_classes + cyclic_eg_classes:
        for d in range(depthSearch):
            QC.update()
    correct_with_mutations = 0
    unknown_acyclic = 0
    unknown_cyclic = 0
    for Q,l in dataset_not_local:
        acyclic_match = False
        Qp = descendFork(Q, maxSteps=depthSearch)
        Qclass = mutationClass(Qp, perms=permutations(Q.n))
        for Ra in acyclic_eg_classes:
            res, Qclass, _ = mutationEquivalent(Qp, Ra.initialQ, depthSearch, Qclass = Qclass, Rclass = Ra, upToIsomorphism=True, returnClasses=True)
            if res is True: 
                acyclic_match = True
                if not l:
                    assert False, f"Quiver {Q.matrix} mutation equivalent to acyclic {Ra.initialQ.matrix} yet is labeled cyclic"
                correct_with_mutations += 1
                break
        if acyclic_match:
            continue
        for Ra in cyclic_eg_classes:
            res, Qclass, _ = mutationEquivalent(Qp, Ra.initialQ, depthSearch, Qclass = Qclass, Rclass = Ra, upToIsomorphism=True, returnClasses=True)
            if res is True:
                if l:
                    assert False, f"Quiver {Q.matrix} mutation equivalent to cyclic {Ra.initialQ.matrix} yet is labeled acyclic"
                correct_with_mutations += 1
                break
        if l:
            unknown_acyclic +=1
        else:
            unknown_cyclic +=1
    c = correct+correct_with_mutations
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f"*Unknown->True: {c+unknown_acyclic} ({100*(c+unknown_acyclic)/D}%)")
    print(f" Unknown->False: {c+unknown_cyclic} ({100*(c+unknown_cyclic)/D}%)")
    
    if not useExhaustive: #do not do the following if exhaustive search is disabled.
        return 

    # baseline 4: Exhaustively search acyclics.
    print()
    print("Comparing to the invariants of all connected acyclics with sufficiently small weights.")
    max_markov = 0
    ns = []
    for Q,l in dataset_not_local:
        m = Q.markov()
        if m and m > max_markov:
            max_markov = m
        if Q.n not in ns:
            ns.append(Q.n)
    print(f"Max Markov Invariant: {max_markov}.")
    acyclic_invs = []
    for n in ns:
        for R in generateAcyclicsBelowMarkov(max_markov, n):
            d = invariantsDictionary(R, includeSurface=useSurface)
            d["Alexander polynomial"] = R.alexanderPolynomial()
            if d not in acyclic_invs:
                acyclic_invs.append(d)
    
    correct_exhaustive = 0 
    unknown_acyclic = 0
    unknown_cyclic = 0
    for Q,l in dataset_not_local:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in acyclic_invs:
            assert not l, f"Invariants {d} are not acyclic yet label is cyclic"
            correct_exhaustive += 1
        else:
            if l: unknown_acyclic+=1
            else: unknown_cyclic+=1

    c = correct+correct_exhaustive
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f"*Unknown->True: {c+unknown_acyclic} ({100*(c+unknown_acyclic)/D}%)")
    print(f" Unknown->False: {c+unknown_cyclic} ({100*(c+unknown_cyclic)/D}%)")

def runBaselineAcyclicLocal(dataset, useExhaustive=False):
    """
    baselines the list of tuples (quiver, label) in dataset against our local-mutation-acyclicity tests.
    if useExhaustive is True, generate all acyclic quivers which might have the same alexander poly (slow!).
    """
    dataset_copy = list(dataset) #could use itertool.tee(dataset, k) and run in parallel for memory savings
    D = len(dataset_copy)

    correct=0
    U_T_correct=0
    U_F_correct=0
    acyclic_count=0
    cyclic_count=0
    # baseline 1: blind attempts. 
    print("Local tests")
    #  3 ways to round: unknown is incorrect; unknown as True, unknown takes False.
    for Q,l in dataset_copy:
        result = mutationAcyclicLocal(Q, matchAlexander=useExhaustive)
        if result is not Unknown:
            assert result == l, f"Dataset and baseline disagree; {Q.matrix} has label {l} but we find {result}."
            correct += 1
        else:
            if l is True:
                U_T_correct+=1
            if l is False:
                U_F_correct+=1
        if l: 
            acyclic_count += 1
        else: 
            cyclic_count += 1
    print(f"data balance (a vs c): {acyclic_count} vs {cyclic_count}")
    print(f" Unknown->incorrect: {correct} ({100*correct/D}%)")
    print(f"*Unknown->True: {correct+U_T_correct} ({100*(correct+U_T_correct)/D}%)")
    print(f" Unknown->False: {correct+U_F_correct} ({100*(correct+U_F_correct)/D}%)")
    
def runBaselineAcyclicWithTraining(dataset, useLocal=True, useSurface=True, train=None, seed=1052026):
    """
    run baseline against dataset using optional training data.
    If no training data is specified, randomly selects ~5% from the dataset using seed.
    May perform quite poorly if dataset is very small or has much more than 20 classes.
    useLocal = True will first use the local mutation acyclicity test. Set to False to disable.
    """

    dataset_copy = list(dataset) #could use itertool.tee(dataset, k) and run in parallel for memory savings
    D = len(dataset_copy)

    # baseline 2: with training data.
    print("Comparing against a train set.")
    if train is None:
        # select a random selection of quivers of each label
        #  we choose to take samples of ~1/20 of the whole dataset; a bit weak if dataset contains few mutations from many classes
        #  provide a comparable train set to your model otherwise.

        acyclic_count=len([d[0] for d in dataset_copy if d[1]])
        cyclic_count=len([d[0] for d in dataset_copy if not d[1]])
        acyclic_egs = random.choices([d[0] for d in dataset_copy if d[1]], k=round(acyclic_count/20))
        cyclic_egs = random.choices([d[0] for d in dataset_copy if not d[1]], k=round(cyclic_count/20))
    else:
        acyclic_egs = [d[0] for d in train if d[1]]
        cyclic_egs = [d[0] for d in train if not d[1]]
    # build useful invariant vectors for each class
    acyclic_invariants = []
    for Q in acyclic_egs:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in acyclic_invariants:
            acyclic_invariants.append(d)

    cyclic_invariants = []
    for Q in cyclic_egs:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in cyclic_invariants:
            cyclic_invariants.append(d)

    print(f"Found {len(cyclic_invariants)} distinct cyclic invariants and {len(acyclic_invariants)} distinct acyclic invariants.")
    overlap = len([d for d in acyclic_invariants if d in cyclic_invariants])
    if overlap > 0:
        print(f"They intersect for {overlap} values.")
    else:
        print("The sets of invariants are mutually distinct.")
    print()

    #correct from local tests (stays 0 if not useLocal)
    correct=0
    #attempt to label the unlabeled with the invariants
    correct_inv_with_cyclic_Alexander = 0
    correct_inv_without_cyclic_Alexander = 0
    undetermined_with_cyclic_Alexander = [0,0] #pair, [#correct label is False, #correct label is True]
    undetermined_without_cyclic_Alexander = [0,0]
    for Q,l in dataset_copy:
        if useLocal:
            r = mutationAcyclicLocal(Q)
            if r is not Unknown:
                assert r == l, f"Dataset and baseline disagree on {Q.matrix}."
                correct += 1
                continue
        d = invariantsDictionary(Q, includeSurface=useSurface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        #match if invariants are in acyclic, in cyclic (without alexander considered), and in cyclic with alexander considered
        match (d in acyclic_invariants, any([all([d[k] == inv[k] for k in d.keys()]) for inv in cyclic_invariants]), d in cyclic_invariants):
            case (True, True, True):
                #ambiguous still
                undetermined_with_cyclic_Alexander[l] += 1 #hack; True=1 as an int while False=0, so this increments the count according to label.
                undetermined_without_cyclic_Alexander[l] += 1
            case (True, False, False):
                if l==True:
                    correct_inv_without_cyclic_Alexander += 1
                    correct_inv_with_cyclic_Alexander += 1
            case (True, True, False):
                undetermined_without_cyclic_Alexander[l] += 1
                if l==True:
                    correct_inv_with_cyclic_Alexander += 1
            case (False, True, True):
                if l==False:
                    correct_inv_without_cyclic_Alexander += 1
                    correct_inv_with_cyclic_Alexander += 1
            case (False, True, False):
                undetermined_with_cyclic_Alexander[l] += 1
                if l==False:
                    correct_inv_without_cyclic_Alexander += 1
            case (False, False, False):
                undetermined_with_cyclic_Alexander[l] += 1
                undetermined_without_cyclic_Alexander[l] += 1
            case (_b, False, True):
                assert False, "Implementation error in heuristic."

    print("Ignoring Alexander Polynomial for cyclic examples:")
    c = correct + correct_inv_without_cyclic_Alexander
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f"*Unknown->True: {c+undetermined_without_cyclic_Alexander[True]} ({100*(c+undetermined_without_cyclic_Alexander[True])/D}%)")
    print(f" Unknown->False: {c+undetermined_without_cyclic_Alexander[False]} ({100*(c+undetermined_without_cyclic_Alexander[False])/D}%)")
    
    print("Using Alexander Polynomials for the cyclic examples:")
    c = correct + correct_inv_with_cyclic_Alexander
    print(f" Unknown->incorrect: {c} ({100*(c)/D}%)")
    print(f" Unknown->True: {c+undetermined_with_cyclic_Alexander[True]} ({100*(c+undetermined_with_cyclic_Alexander[True])/D}%)")
    print(f" Unknown->False: {c+undetermined_with_cyclic_Alexander[False]} ({100*(c+undetermined_with_cyclic_Alexander[False])/D}%)")

def runBaselineMatchSets(dataset, classReps, useSurface=True):
    """ 
    For each quiver in classReps, create a set of invariants and match the quivers in dataset accordingly.
    dataset has tuples of quivers and labels, indicating the index of the class the quiver is from.
    runs both with and without the alexander polynomial for (possibly mutation cyclic) quivers.
    The alexander polynomial of an acyclic or 3-vert classRep is used for both cases.
    """
    D=len(dataset)
    invariantsWithAlex = []
    invariantsWithoutAlex = []
    totallyProperRep = []
    for Q in classReps:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        da = copy.deepcopy(d)
        da["Alexander polynomial"] = Q.alexanderPolynomial()
        if Q.acyclic() is True or Q.n <= 3:
            totallyProperRep.append(True) 
        else:
            totallyProperRep.append(False)
        invariantsWithoutAlex.append(d)
        invariantsWithAlex.append(da)

    correctWithout = 0
    correctWith = 0
    guessWithout = 0.0
    guessWith = 0.0
    for Q,l in dataset:
        d = invariantsDictionary(Q, includeSurface=useSurface)
        da = copy.deepcopy(d)
        da["Alexander polynomial"] = Q.alexanderPolynomial()
        optionsWithout = []
        optionsWith = []
        for i in range(len(classReps)):
            if invariantsWithoutAlex[i]!=d:
                continue
            if totallyProperRep[i]:
                if invariantsWithAlex[i]!=da:
                    continue
            optionsWithout.append(i)
            if invariantsWithAlex[i]!=da:
                continue
            optionsWith.append(i)
        assert l in optionsWithout, f"Dataset and baseline disagree on {Q.matrix}, expected {l} in {optionsWithout}"
        if len(optionsWithout)==1:
            correctWith+=1
            correctWithout+=1
            continue
        else:
            guessWithout += 1/len(optionsWithout)
        if len(optionsWith)==1:
            correctWith+=1
            continue
        if len(optionsWith)==0:
            # reaching this branch implies we had multiple to guess with optionsWith; use the same guess.
            guessWith += 1/len(optionsWithout)
        else:
            guessWith += 1/len(optionsWith)

    print(f"Ignoring Alexander Polynomial for non-totally-proper reps: {correctWithout} ({100*(correctWithout)/D}%)")
    print(f"... and uniform guessing: {correctWithout + guessWithout} ({100*(correctWithout + guessWithout)/D}%)")
    
    print(f"Using Alexander Polynomials for the cyclic examples: {correctWith} ({100*(correctWith)/D}%)")
    print(f"... and uniform guessing: {correctWith + guessWith} ({100*(correctWith + guessWith)/D}%)")
    


