from Quiver import *

# Contains algorithms and helper functions mirroring how we would solve 
# mutation equivalence or mutation acyclicity
#  includes variations with more/less compute, certainty, and data. 

def invariants_dict(Q, include_surface=True):
    # return a dict of mutation invariants of Q. 
    invs = dict()
    invs["rank"] = Q.n
    invs["connected"] = Q.connected()
    invs["B rank"] = Q.B_rank()
    invs["Casal's Det"] = Q.casals_det()
    invs["Seven's Congruence"] = Q.seven_congruence()
    invs["gcd"] = Q.gcd_vector() 
    if include_surface: invs["surf quiver"] = surface_quiver(Q)[0] 
    #note, block decomp is not an invariant but being surface type is
    return invs

def mutation_equivalent_local(Q,R, up_to_ismorphism=False):
    # return what we can determine about the mutation-equivalence of Q and R
    # Does not mutate either quiver. 
    #  if up_to_ismorphism is True, compares Q against all quivers isomorphic to R
    if up_to_ismorphism:
        all_known = True
        for p in permutations(R.n()):
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
        if Q.alexander_poly() != R.alexander_poly():
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
    if Q.hasMutCyclicSubquiver():
        return False
    if Q.markov() is False: #includes Q contains a vortex
        return False
    if Q.markov() < 0 or (Q.connected() and Q.markov() < Q.n-1):
        return False
    if Q.determinant()**2 > (2*Q.markov())**Q.n:
        return False
    if match_alexander:
        for R in generate_acyclics_from_Alexander(Q.alexander_poly()):
            if mutation_equivalent_local(Q,R) is not False:
                return Unknown #Heuristically, likely to be True
        return False
    else:
        return Unknown


def run_benchmark_acyclic(dataset, train=None, use_surface=True, use_exhaustive=False, seed=2842026):
    # Runs several variants of our benchmark against the given dataset, and prints the accuracies achieved.
    #  dataset: iterable of tuples, (quiver, correct_label), quiver is a quiver object and correct_label is a bool
    #  train (optional): iterable of (quiver, correct_label) pairs 
    #   (train may be used to guess other labels, if None then some benchmarks will sample examples from dataset)
    #  use_surface: if True, computes the block decomposition for invariants.
    #  use_exhaustive: if True, enumerates all acyclic quivers with relevant size.
    #  seed: used for random in selecting subset of data, irrelevant if train is not None.

    dataset_copy = list(dataset) #could use itertool.tee(dataset, k) and run in parallel for memory savings
    D = len(dataset_copy)
    dataset_unknowns = []
    acyclic_count = 0
    cyclic_count = 0
    random.seed(seed)
    # benchmark 1: blind attempts. 
    #  3 ways to round: unknown is incorrect; unknown as True, unknown takes False.
    for Q,l in dataset_copy:
        result = mutation_acyclic_local(Q, match_alexander==use_exhaustive)
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

    print("Local tests")
    print(f" Unknown->incorrect: {correct} ({100*correct/D}%)")
    print(f"*Unknown->True: {correct+U_T_correct} ({100*(correct+U_T_correct)/D}%)")
    print(f" Unknown->False: {correct+U_F_correct} ({100*(correct+U_F_correct)/D}%)")

    # benchmark 2: with training data.
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
    def get_useful_invariants(egs):
        # return a list of distinct invariant dicts. 
        #  for each quiver in egs, compute the invariant dict and add it to the result if it is not implied by an existing vector.
        #  remove any dicts it implies
        
        #def implies(T1, T2):
        #    # generous ternary implies; meaning T1 has more info than T2.
        #    if T2 is Unknown:
        #        return True
        #    if T1 is Unknown:
        #        return False
        #    return T1 == T2

        
        result_dicts = []
        for Q in egs:
            d = invariants_dict(Q, include_surface=use_surface)
            d["Alexander poly"] = Q.alexander_poly()
            if d not in result_dicts:
                result_dicts.append(d)
        return result_dicts
    acyclic_invariants = []
    cyclic_invariants = []



