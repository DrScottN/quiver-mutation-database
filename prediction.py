from Quiver import *
import random

# Contains algorithms and helper functions mirroring how we would solve 
# mutation equivalence or mutation acyclicity
#  includes variations with more/less compute, certainty, and data. 

def invariants_dict(Q, include_surface=True):
    # return a dict of mutation invariants of Q. 
    invs = dict()
    invs["rank"] = Q.n
    invs["connected"] = Q.connected()
    invs["B rank"] = matrix_rank(Q.matrix)
    invs["Casal's Det"] = Q.casalsDeterminant()
    invs["Seven's Congruence"] = Q.seven_congruence()
    invs["gcd"] = Q.gcdVector() 
    if include_surface: invs["surf quiver"] = bool(surface_quiver(Q))
    #note, block decomp is not an invariant but being surface type is
    return invs

def mutation_equivalent_local(Q,R, up_to_isomorphism=False):
    # return what we can determine about the mutation-equivalence of Q and R
    # Does not mutate either quiver. 
    #  if up_to_isomorphism is True, compares Q against all quivers isomorphic to R
    if up_to_isomorphism:
        all_known = True
        for p in permutations(R.n):
            result = mutation_equivalent_local(Q, isomorphicQuiver(R, p))
            if result is True:
                return result
            if result is Unknown:
                all_known = False
        if all_known:
            return False
        return Unknown
    
    if invariants_dict(Q) != invariants_dict(R):
        return False
    if Q.acyclic() or R.acyclic():
        if Q.alexanderPolynomial() != R.alexanderPolynomial():
            return False
    return Unknown


def mutation_acyclic_local(Q, match_alexander=False):
    # return what we can determine about the mutation-acyclicity of Q
    # Does not mutate Q, move to forkless part first for best results.
    #  if match_alexander=True, then generates all acyclic quivers with the same alexander poly
    #  and checks if they could be mutation equivalent (based on local data).
    #   This can be extremely slow for large quivers.
    # suggested convention: round Unknown to match_alexander
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
    if match_alexander:
        for R in generate_acyclics_from_Alexander(Q.alexanderPolynomial()):
            if mutation_equivalent_local(Q,R) is not False:
                return Unknown #Heuristically, likely to be True
        return False
    else:
        return Unknown

def mutation_equivalent(Q, R, distance, up_to_isomorphism=False, return_classes=False, Qclass=None, Rclass=None):
    """Determines if Q is isomorphic to R using ~distance mutations, along with the invariants.
    returns True if they are equivalent, False if they are not, and Unknown if we cannot determine either way.
    
    up_to_isomorphism=True considers the problem up to isomorphism
    return_classes=True returns a tuple (Result, Mutation_class_of_Q, Mutation_class_of_R)
    Q_class, R_class may be provided for efficiency; distance is ignored if a class is given.
    attempts to replace Q,R with non-forks in distance mutations if class is not provided. 
     This may result in more than distance mutations being applied."""
    if up_to_isomorphism:
        permsQ = permutations(Q.n)
        permsR = permutations(R.n)
    else:
        permsQ = [list(range(Q.n))]
        permsR = [list(range(R.n))]

    update = False
    if Qclass is None:
        Q = descend_fork(Q, maxSteps=distance)
        Qclass = mutationClass(Q, perms=permsQ)
        update = True
    if Rclass is None:
        R = descend_fork(R, maxSteps=distance)
        Rclass = mutationClass(R, perms=permsR)
        update = True

    local_check = mutation_equivalent_local(Q, R, up_to_isomorphism=up_to_isomorphism)
    if local_check is not Unknown:
        return (local_check, Qclass, Rclass) if return_classes else local_check
    
    #we may now assume same rank, etc.

    if update:
        distance = (distance+1)//2
        for d in range(distance):
            R_Class.update()
            Q_class.update()
    
    #check if overlap found
    if union := Rclass.union(Qclass):
        return (True, union, union) if return_classes else True

    #check if invariants are consistent
    if Rclass.agrees(Qclass) is False: #not implemented
        return (False, Qclass, Rclass) if return_classes else False

    return (Unknown, Qclass, Rclass) if return_classes else Unknown




def runBenchmarkAcyclic(dataset, depthSearch=1, train=None, use_surface=True, use_exhaustive=False, seed=2842026):
    """Runs several variants of our benchmark against the given dataset, and prints the accuracies achieved.
        dataset: iterable of tuples, (quiver, correct_label), quiver is a quiver object and correct_label is a bool
        train (optional): iterable of (quiver, correct_label) pairs 
        (train may be used to guess other labels, if None then some benchmarks will sample examples from dataset)
        use_surface: if True, computes the block decomposition for invariants.
        use_exhaustive: if True, enumerates all acyclic quivers with relevant size.
        seed: used for random in selecting subset of data, irrelevant if train is not None."""

    dataset_copy = list(dataset) #could use itertool.tee(dataset, k) and run in parallel for memory savings
    D = len(dataset_copy)
    dataset_not_local = []
    acyclic_count = 0
    cyclic_count = 0
    random.seed(seed)

    correct=0
    U_T_correct=0
    U_F_correct=0
    # benchmark 1: blind attempts. 
    print("Local tests")
    #  3 ways to round: unknown is incorrect; unknown as True, unknown takes False.
    for Q,l in dataset_copy:
        result = mutation_acyclic_local(Q, match_alexander=use_exhaustive)
        if result is not Unknown:
            assert result == l, f"Dataset and benchmark disagree; {Q.matrix} has label {l} but we find {result}."
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
    
    # benchmark 2: with training data.
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
        d = invariants_dict(Q, include_surface=use_surface)
        d["Alexander polynomial"] = Q.alexanderPolynomial()
        if d not in acyclic_invariants:
            acyclic_invariants.append(d)

    cyclic_invariants = []
    for Q in cyclic_egs:
        d = invariants_dict(Q, include_surface=use_surface)
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
        d = invariants_dict(Q, include_surface=use_surface)
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

    # benchmark 3: Apply mutations
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
        Qp = descend_fork(Q, maxSteps=depthSearch)
        Qclass = mutationClass(Qp, perms=permutations(Q.n))
        for Ra in acyclic_eg_classes:
            res, Qclass, _ = mutation_equivalent(Qp, Ra.initialQ, depthSearch, Qclass = Qclass, Rclass = Ra, up_to_isomorphism=True, return_classes=True)
            if res is True: 
                acyclic_match = True
                if not l:
                    assert False, f"Quiver {Q.matrix} mutation equivalent to acyclic {Ra.initialQ.matrix} yet is labeled cyclic"
                correct_with_mutations += 1
                break
        if acyclic_match:
            continue
        for Ra in cyclic_eg_classes:
            res, Qclass, _ = mutation_equivalent(Qp, Ra.initialQ, depthSearch, Qclass = Qclass, Rclass = Ra, up_to_isomorphism=True, return_classes=True)
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
    
    if not use_exhaustive: #do not do the following if exhaustive search is disabled.
        return 

    # benchmark 4: Exhaustively search acyclics.
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
            d = invariants_dict(R, include_surface=use_surface)
            d["Alexander polynomial"] = R.alexanderPolynomial()
            if d not in acyclic_invs:
                acyclic_invs.append(d)
    
    correct_exhaustive = 0 
    unknown_acyclic = 0
    unknown_cyclic = 0
    for Q,l in dataset_not_local:
        d = invariants_dict(Q, include_surface=use_surface)
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
    



