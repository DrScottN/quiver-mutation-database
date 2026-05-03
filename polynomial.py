
class polynomial():
    # need to make these deal with multiple variables
    # as such, coefficients will now be a list of tuples, where the first entry is the coefficient and the remaining are the powers
    # where the order of the tuples is given by the variable order
    # need to let coefficients be a stream eventually to deal with formal power series
    def __init__(self, coefficients, variables = ('x',)) -> None:
        enum_vars = sorted(enumerate(variables), key=lambda x: x[1])
        self.variables = [x[1] for x in enum_vars]
        self.numVars = len(variables)
        self.simplified = False

        if isinstance(coefficients, list):
            if len(coefficients)>0 and not isinstance(coefficients[0], tuple):
                coefficients = [(c, i) for i, c in enumerate(coefficients) if c!=0]
            else:
                coefficients = [(c[0],) + tuple(c[1+x[0]] for x in enum_vars) for c in coefficients if c[0]!=0]
            self.coefficients = coefficients
            self.__lex__() # Sorts them in the lexigraphical order and adds terms of the same monomial
            self.totalDegree = max(sum(c[1:]) for c in self.coefficients+[[0,0]])
            if len(self.variables) == 1:
                self.degree = self.totalDegree
            
            
            self.localDegrees = {v : max(c[i+1] for c in self.coefficients+[[0]+[0]*len(self.variables)]) for i, v in enumerate(self.variables)}
        else:
            raise Exception("Not a valid input")

    def __str__(self) -> str:
        s = " + ".join([f"({c[0]}){''.join(['^'.join([self.variables[i], str(j)]) for i,j in enumerate(c[1:]) if j != 0])}" for c in self.coefficients if c[0] != 0])
        return s

    def __varseq__(self, other):
        # Returns true if two polynomials have the same variables
        a = set(self.variables)
        b = set(other.variables)

        return a == b

    def __varchange__(self, newVars):
        # Modifies the coefficients to match the new variables and returns a new polynomial
        # newVars is a tuple as well
        variables = {v : i+1 for i, v in enumerate(self.variables)}
        totalVars = tuple(set(self.variables) | set(newVars))
        coeff = [(c[0],) + tuple([c[variables[v]] if v in variables else 0 for v in totalVars]) for c in self.coefficients]

        return polynomial(coeff, totalVars)

    def __lex__(self, reverse = True):
        # reformats the polynomial to be in reverse lexigraphical order
        # also builds a dictionary of coefficients, which is sometimes easier
        self.coeffDict = {}
        for c in self.coefficients:
            if c[0] != 0:
                self.coeffDict[c[1:]] = c[0]

        #zero = (0,) * (len(self.coefficients[0])-1)

        #if zero not in self.coeffDict:
        #    self.coeffDict[zero] = 0

        self.coeffDict = dict(sorted(self.coeffDict.items())) # Should automatically sort
        newCoeff = [(v,) + k for k,v in self.coeffDict.items()]
        self.coefficients = newCoeff if not reverse else newCoeff[::-1]   

    def __eq__(self, other):
        if isinstance(other, int): #or isinstance(other, Rational):
            return self.totalDegree == 0 and ((len(self.coefficients)==0 and other==0) or self.coefficients[0][0] == other)
        elif isinstance(other, polynomial):
            if not self.__varseq__(other):
                newVars = set(self.variables + other.variables)
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
        newVars = self.variables
        if isinstance(other, int): #or isinstance(other, Rational):
            deg_zero = (0,)*len(self.variables)
            coefficients = [c for c in self.coefficients]
            found=False
            for i,c in enumerate(self.coefficients):
                if c[1:] == deg_zero:
                    coefficients[i] = (c[0]+other,) + deg_zero
                    found=True
                    break
            if not found:
                coefficients = [(other,) + deg_zero] + coefficients
        elif isinstance(other, polynomial):
            if not self.__varseq__(other):
                newVars = tuple(sorted(set(other.variables + newVars)))
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

        return polynomial(coefficients, variables=newVars)
    
    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        coefficients = [(-c[0],) + c[1:] for c in self.coefficients]

        return polynomial(coefficients, self.variables)
    
    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if self == 0 or other == 0:
            return polynomial([0], variables=self.variables)
        elif not isinstance(other, polynomial):
            coefficients = [(other * c[0],) + c[1:] for c in self.coefficients]
            newVars = self.variables
        else:
            coeffDict = dict()
            if not self.__varseq__(other):
                newVars = tuple(sorted(set(self.variables + other.variables)))
                p = self.__varchange__(newVars)
                q = other.__varchange__(newVars)
            else:
                p = self
                q = other
                newVars=self.variables

            for c in p.coefficients:
                for d in q.coefficients:
                    monomial = tuple([c[i] + d[i] for i in range(1,len(c))])
                    if monomial not in coeffDict:
                        coeffDict[monomial] = 0

                    coeffDict[monomial] += c[0] * d[0]
                    
            coefficients = [(v,) + k for k,v in coeffDict.items()]

        r = polynomial(coefficients, variables = newVars)

        return r
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if not isinstance(other, int):
            Exception("Can only exponentiate a polynomial with an integer")
        elif other < 0:
            Exception("Can only exponentiate a polynomial with a positive integer")
        elif other == 0:
            return 1
        elif other == 1:
            return self

        a = other // 2
        b = other - a # will be >= a
        #print(self, a, b)
        return pow(self, a)*pow(self, b) # This should be a faster exponentiation

    def __truediv__(self, other):
        """ Exact division of self by other """
        if other == 0:
            Exception("Division by 0")
        if self.numVars > 1 or (isinstance(other, polynomial) and other.numVars > 1):
            Exception("Unsupported: divide with multivariate polynomial.")
        if self==0:
            return 0
        if isinstance(other ,int):
            otherLead = (other, 0)
        else:
            otherLead = [x for x in other.coefficients if x[1]==other.totalDegree][0]
        selfLead = [x for x in self.coefficients if x[1]==self.totalDegree][0]
        if otherLead[1] > selfLead[1]:
            raise Exception("not exact division")
        if selfLead[0] % otherLead[0] != 0:
            raise Exception("Not exact division")
        div = (selfLead[0]//otherLead[0])*(polynomial([0,1], variables=self.variables)**(self.totalDegree - otherLead[1]))
        #print(self, " ", other, "div:", div)
        return div + (self - div*other)/other

    def coefficientList(self):
        """ return list of coefficients for 1 variable poly.
        self is polynomial(self.coefficientList()) """
        if self.numVars!=1: raise Exception("coefficient List supported only for one variable polys")
        res = [0]*self.totalDegree
        for c in self.coefficients:
            res[c[1]] = c[0]
        return res
           
    def eval(self, values):
        # Evalutes the polynomial at variables = values
        s = 0
        
        for c in self.coefficients:
            #print(f"{s}, {c * (val ** i)}")
            s += c[0] * sum(values[i] ** j for i, j in enumerate(c[1:]))

        return s
