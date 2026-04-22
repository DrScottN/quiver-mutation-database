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
    #note, block decomp is not an invariant but the surface type is.
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