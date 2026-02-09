# File to handle all the code for producing labels from a given quiver
import Quiver

# List of all possible labels, copied from our overleaf file
labels = ["mutation class",                     # MutationClass(Quiver) 
          "mutation class up to isomorphism",   # [MutationClass(iso(Quiver)) for iso in isomorphisms]
          "breadth",                            # Integer
          "mutation-acyclic",                   # Boolean or None for unknown 
          "surface quiver",                     # Boolean or None for unknown
          "plabic quiver",                      # Boolean or None for unknown
          "finite type",                        # Boolean or None for unknown
          "mutation-finite",                    # Boolean or None for unknown 
          "finite forkless part",               # Boolean or None for unknown
          "finite pre-forkless part",           # Boolean or None for unknown
          "mutation-abundant",                  # Boolean or None for unknown
          "locally acyclic",                    # Boolean or None for unknown
          "mutation-complete",                  # Boolean or None for unknown
          "mutation vortex-free",               # Boolean or None for unknown
          "reddening sequence",                 # Boolean or None for unknown
          "greenening sequence",                # Boolean or None for unknown
          "totally proper",                     # Boolean or None for unknown
          "contains vortex",                    # Boolean
          "connected",                          # Boolean
          "weight gcd vectors",                 # Integer arrary
          "matrix rank",                        # Integer
          "quiver size",                        # Integer
          "congruence invariant"]               # Integer
# Will need to add more, but this should be a good start for figuring out labels

# Functions for checking the labels

def isAcyclic(Q):
    return Q.acyclic()

def isMutationAcyclic(mutClass):
    return mutClass.acyclic

def isMutationFinite(mutClass):
    return mutClass.finite


# Label function

def label(Q):
    # Returns a dictionary of labels associated to the given quiver Q
    # Not sure how to deal with mutation class yet

    labelled = {l : None for l in labels}

    if Q.acyclic():
        labelled["mutation-acyclic"] = True

    return labelled

if __name__ == "__main__" :
    Q = Quiver.boxQuiver(2,2)
    print(label(Q))


