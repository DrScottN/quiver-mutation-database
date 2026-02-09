
class polynomial():
    # need to make these deal with multiple variables
    # as such, coefficients will now be a list of tuples, where the first entry is the coefficient and the remaining are the powers
    # where the order of the tuples is given by the variable order
    # need to let coefficients be a stream eventually to deal with formal power series
    def __init__(self, coefficients, vars = ('x',)) -> None:
        self.vars = vars
        self.numVars = len(vars)
        self.simplified = False

        if isinstance(coefficients, list):
            if not isinstance(coefficients[0], tuple):
                coefficients = [(c, i) for i, c in enumerate(coefficients)]

            self.coefficients = coefficients
            self.__lex__() # Sorts them in the lexigraphical order and adds terms of the same monomial
            self.totalDegree = max(sum(c[1:]) for c in self.coefficients)
            if len(coefficients[0][1:]) == 1:
                self.degree = self.totalDegree
            
            
            self.localDegrees = {v : max(c[i+1] for c in self.coefficients) for i, v in enumerate(self.vars)}
        else:
            raise Exception("Not a valid input")

    def __str__(self) -> str:
        s = " + ".join([f"({c[0]}){''.join(['^'.join([self.vars[i], str(j)]) for i,j in enumerate(c[1:]) if j != 0])}" for c in self.coefficients if c[0] != 0])
        return s

    def __varseq__(self, other):
        # Returns true if two polynomials have the same variables
        a = set(self.vars)
        b = set(other.vars)

        return a == b

    def __varchange__(self, newVars):
        # Modifies the coefficients to match the new variables and returns a new polynomial
        # newVars is a tuple as well
        vars = {v : i+1 for i, v in enumerate(self.vars)}
        totalVars = list(set(self.vars) | set(newVars))
        coeff = [(c[0],) + tuple([c[vars[v]] if v in vars else 0 for v in totalVars]) for c in self.coefficients]

        return polynomial(coeff, totalVars)

    def __lex__(self, reverse = True):
        # reformats the polynomial to be in reverse lexigraphical order
        # also builds a dictionary of coefficients, which is sometimes easier
        self.coeffDict = {}
        for c in self.coefficients:
            if c[0] != 0:
                self.coeffDict[c[1:]] = c[0]

        zero = (0,) * (len(self.coefficients[0])-1)

        if zero not in self.coeffDict:
            self.coeffDict[zero] = 0

        self.coeffDict = dict(sorted(self.coeffDict.items())) # Should automatically sort
        newCoeff = [(v,) + k for k,v in self.coeffDict.items()]
        self.coefficients = newCoeff if not reverse else newCoeff[::-1]   

    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, Rational):
            return self.totalDegree == 0 and self.coefficients[0][0] == other
        elif isinstance(other, polynomial):
            if not self.__varseq__(other):
                newVars = set(self.vars + other.vars)
                p = self.__varchange__(newVars)
                q = other.__varchange__(newVars)
            else:
                p = self
                q = other

            if p.totalDegree != q.totalDegree:
                return False
            
            for c, k in zip(p.coefficients, q.coefficients):
                if c != k:
                    return False
                
            return True
        else:
            raise Exception("What are you comparing?")         

    def __add__(self, other):
        newVars = self.vars
        if isinstance(other, int) or isinstance(other, Rational):
            coefficients = [(c[0] + other,) + c[1:] for c in self.coefficients]
        elif isinstance(other, polynomial):
            if not self.__varseq__(other):
                newVars = list(set(self.vars + newVars))
                p = self.__varchange__(newVars)
                q = other.__varchange__(newVars)
            else:
                p = self
                q = other

            pDict = p.coeffDict
            qDict = q.coeffDict
            newCoeff = []

            for k, v in pDict.items():
                added = v
                if k in qDict:
                    added += qDict[k]
                newCoeff.append((added,) + k)

            for k, v in qDict.items():
                if k in pDict:
                    continue
                newCoeff.append((v,) + k)

            coefficients = newCoeff[:]
        else:
            raise Exception("Unsupported addition")

        return polynomial(coefficients, vars=newVars)
    
    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        coefficients = [(-c[0],) + c[1:] for c in self.coefficients]

        return polynomial(coefficients, self.vars)
    
    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if self == 0 or other == 0:
            return polynomial([0], var=self.vars)
        elif not isinstance(other, polynomial):
            coefficients = [(other * c[0],) + c[1:] for c in self.coefficients]
            newVars = self.vars
        else:
            coeffDict = dict()
            if not self.__varseq__(other):
                newVars = list(set(self.vars + other.vars))
                p = self.__varchange__(newVars)
                q = other.__varchange__(newVars)
            else:
                p = self
                q = other

            for c in p.coefficients:
                for d in q.coefficients:
                    monomial = tuple([c[i] + d[i] for i in range(1,len(c))])
                    if monomial not in coeffDict:
                        coeffDict[monomial] = 0

                    coeffDict[monomial] += c[0] * d[0]
                    
            coefficients = [(v,) + k for k,v in coeffDict.items()]

        r = polynomial(coefficients, vars = newVars)

        return r
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def eval(self, values):
        # Evalutes the polynomial at vars = values
        s = 0
        
        for c in self.coefficients:
            #print(f"{s}, {c * (val ** i)}")
            s += c[0] * sum(values[i] ** j for i, j in enumerate(c[1:]))

        return s
