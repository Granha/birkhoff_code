import sys
import fractions
import decimal

#################################
#  Decimal conversion function  #
#################################

def toDecimal(number, precision=21, *, rounding=decimal.ROUND_CEILING, Emin=decimal.MIN_EMIN, Emax=decimal.MAX_EMAX):
    with decimal.localcontext(decimal.Context(prec=precision, rounding=rounding, Emin=Emin, Emax=Emax)):
        if isinstance(number, fractions.Fraction):
            return decimal.Decimal(number.numerator) / decimal.Decimal(number.denominator)
        return +decimal.Decimal(number) # The unary plus is used to convert to the precision passed

################################################
#  Callback functions and function generators  #
################################################

def dummy(**kwargs):
    """
    Does not do anything.
    """
    pass

def joinFunctions(*args):
    """
    Receives a list of functions and returns a function that executes them in order on the same **kwargs.
    The sequence is halted at the first function that returns a not None value and such value
    is returned; otherwise returns None.
    """
    def joined(**kwargs):
        for f in args:
            r = f(**kwargs)
            if r != None:
                return r
        return None
    return joined

def testItMod(mod=1):
    """
    Returns a function that tests whether **kwargs has 'it' with a value divisible by mod or if it has 'solved'
    with value True.
    The returned function returns None if either this is True or 'it' is None or does not exist
    in **kwargs; otherwise returns False.
    """
    def function(**kwargs):
        it = kwargs.get('it', None)
        solved = kwargs.get('solved', None)
        if it == None or it % mod == 0 or solved==True:
            return None
        return False
    return function

def printer(name, formatstr='{0} =\n{1}', noneFormat='{0} = None'):
    """
    Returns a function that prints element name from **kwargs using format formatstr and using
    format noneFormat if the value is None.
    If **kwargs does not contain name, treats as if it had value None.
    """
    def function(**kwargs):
        item = kwargs.get(name, None)
        if item == None:
            print(noneFormat.format(name, item), end='')
        else:
            print(formatstr.format(name, item), end='')
    return function

def dictprinter(name, strdict):
    """
    Returns a function that prints strdict[value] if name is in **kwargs with value value.
    If strdict[value] is None or strdict does not contain value, does not print anything.
    If **kwargs does not contain name, treats as if it had value None.
    """
    def function(**kwargs):
        item = kwargs.get(name, None)
        string = strdict.get(item, None)
        if string != None:
            print(string, end='')
    return function


def c0printer(precision=21):
    """
    Returns a function that prints element 'c0' of **kwargs as a Decimal with precision precision, if it exists.
    If 'maximization' is in **kwargs and is False and 'phase' in **kwargs is not 'Feasibility', prints -c0
    instead.
    """
    def function(**kwargs):
        c0 = kwargs.get('c0', None)
        if c0 != None:
            if kwargs.get('maximization', True) == False and kwargs.get('phase', None) != 'Feasibility':
                c0 = -c0
            print('%s ' % str(toDecimal(c0, precision)), end='')
    return function

def c0diffPrinter(precision=3, formatstr='(step: %s) '):
    """
    Returns a function that prints the difference between element 'c0' of **kwargs using format formatstr and
    its last value as a Decimal with precision precision, if it exists. Does not print anything if 'it' in
    **kwargs is 0 (but still reads the value of c0 to be used for next difference).
    """
    def function(**kwargs):
        c0 = kwargs.get('c0', None)
        if c0 != None:
            if kwargs.get('it', None) != 0:
                print(formatstr % str(toDecimal(c0 - function.lastc0, precision)), end='')
            function.lastc0 = c0

    return function
    

def generateModPrinter(mod=1, precision=21, diffprecision=3):
    """
    Returns a printer of every mod iterations that prints phase, repetition, iteration and objective as Decimal
    with precision precision.
    """
    return joinFunctions(testItMod(mod),
                         printer('phase', '{1} Phase, ', ''),
                         printer('rep', 'Repetition {1}: ', ''),
                         printer('it', 'Iteration {1}: ', ''),
                         printer('c0', 'Objective: ', ''),
                         c0printer(precision),
                         c0diffPrinter(diffprecision),
                         dictprinter('solved', {True:'Solved! '}),
                         dictprinter('optimal', {True:'Optimal point found! '}),

                         lambda **kwargs : print() #end of line
                         )

def generateStateSaver(d):
    """
    Returns a function that saves **kwargs in d. (Old elements of d are deleted and all saves are shallow.)
    """
    def function(**kwargs):
        d.clear()
        for (i,v) in kwargs.items():
            d[i] = v
    return function

def printAll(**kwargs):
    """
    Prints **kwargs.
    """
    print(kwargs)

def sizeOfFraction(a):
    """
    Returns an integer that approximately captures the size of a in memory (in bytes?)
    """ 
    return sys.getsizeof(a.numerator) + sys.getsizeof(a.denominator)

def _traverseNested(obj, *, top=False):
    """
    Returns a generator to pairs (keylist, val) that completely traverses obj and its nested lists and dicts.
    Each (keylist, val) is of the form keylist = [i1, ..., it] such that
    - if top=False, then obj[it]...[i1] = val.
    - if top=True, then obj[i1]...[it] = val.
    """
    if type(obj) is dict:
        for key in obj:
            for (keylist, val) in _traverseNested(obj[key]):
                keylist.append(key)
                if top:
                    keylist.reverse()
                yield (keylist, val)
    elif type(obj) is list:
        for key in range(len(obj)):
            for (keylist, val) in _traverseNested(obj[key]):
                keylist.append(key)
                if top:
                    keylist.reverse()
                yield (keylist, val)
    else:
        yield ([], obj)

        
def traverseNested(obj):
    """
    Returns a generator to pairs (keylist, val) that completely traverses obj and its nested lists and dicts.
    Each (keylist, val) is of the form keylist = [i1, ..., it] such that obj[i1]...[it] = val.
    """
    return _traverseNested(obj, top=True)

def generateMaxSizeTracker(initialSize=100, *, inspector=sizeOfFraction, items=['A', 'c', 'b', 'c0'],
                           formatstr='Size {0} at {1}.\n'):
    """
    Returns a function that tracks the maximum size occupied by elements in the system. It reports only when
    the maximum size increases. It only starts reporting when the maximum size is larger than initialSize.
    The size of each element is computed using inspector and the items inspected are listed in items.
    The report is formatted according to formatstr. 
    """
    def tracker(**kwargs):
        for item in items:
            obj = kwargs.get(item, None)
            if obj == None:
                continue
            for (keylist, val) in traverseNested(obj):
                size = inspector(val)
                if size > tracker.maxSize:
                    tracker.maxSize = size
                    
                    atstr = str(item)
                    for key in keylist:
                        atstr += '[' + str(key) + ']'

                    print(formatstr.format(size, atstr, val), end='')

    tracker.maxSize = initialSize
    return tracker

def generateMaxSizeTrackerRevealer(initialSize=100, *, inspector=sizeOfFraction, items=['A', 'c', 'b', 'c0']):
    """
    Returns a function that tracks the maximum size occupied by elements in the system. It reports only when
    the maximum size increases. It only starts reporting when the maximum size is larger than initialSize.
    The size of each element is computed using inspector and the items inspected are listed in items.
    It also prints the largest element when the size increases.
    """
    return generateMaxSizeTracker(initialSize, inspector=inspector, items=items,
                                  formatstr='Size {0} at {1}:\nValue = {2}\n')

def generateAverageSizeTracker(*, inspector=sizeOfFraction, items=['A', 'c', 'b', 'c0'], precision=9):
    """
    Returns a function that tracks the average size occupied by elements in the system. 
    The size of each element is computed using inspector and the items inspected are listed in items.
    """
    def tracker(**kwargs):
        total = 0
        number = 0
        for item in items:
            obj = kwargs.get(item, None)
            if obj == None:
                continue
            for (keylist, val) in traverseNested(obj):
                size = inspector(val)
                total += size
                number += 1
        if number > 0:
            print('Average size: %s' % toDecimal(total / number, precision))
    return tracker

def generateTotalSizeTracker(*, inspector=sizeOfFraction, items=['A', 'c', 'b', 'c0']):
    """
    Returns a function that tracks the total size occupied by elements in the system. 
    The size of each element is computed using inspector and the items inspected are listed in items.
    """
    def tracker(**kwargs):
        total = 0
        for item in items:
            obj = kwargs.get(item, None)
            if obj == None:
                continue
            for (keylist, val) in traverseNested(obj):
                size = inspector(val)
                total += size
        print('Total size: %s' % str(total))
    return tracker

def computeDensity(a):
    """
    Returns the density (i.e., fraction of non-zero elements) of a list of numbers or dict whose values
    are numbers.
    """
    if type(a) is list:
        return sum(1 for val in A[i] if val != 0) / len(A[i])
    if type(a) is dict:
        return sum(1 for val in A[i].values() if val != 0) / len(A[i])

def generateDensityTracker(*, precision=9):
    """
    Returns a function that tracks the minimum and average density across the system.
    The precision of the printed density is given by precision.
    """
    def tracker(**kwargs):
        minimum = 2.0
        sparsest = None
        total = 0.0
        number = 0

        A = kwargs.get('A', None)
        if A != None:
            for i in range(len(A)):
                if len(A[i]) > 0:
                    density = computeDensity(A[i])
                    if density < minimum:
                        sparsest = 'A[%d]' % i
                        minimum = density
                    total += density
                    number += 1

        c = kwargs.get('c', None)
        if c != None and len(c) > 0:
            density = computeDensity(c)
            if density < minimum:
                sparsest = 'c'
                minimum = density
            total += density
            number += 1

        b = kwargs.get('b', None)
        if b != None and len(b) > 0:
            density = computeDensity(b)
            if density < minimum:
                sparsest = 'b'
                minimum = density
            total += density
            number += 1
        if number > 0:
            print('Average density: %s' % str(toDecimal(total / number, precision)))
        if minimum != 2.0:
            print('Minimum density: %s at %s' % (str(toDecimal(minimum, precision)), sparsest))
    return tracker


#############################
#  Pivot choice strategies  #
#############################

def basicStrategy(c, A, b, c0, base):
    """
    Simple naive strategy for pivot choice.

    Returns the pivot (i0, pivot). Returns None if the current solution is optimal.
    """
    (cmax, pivot) = max(((c[var], var) for var in c))

    if cmax <= 0:
        return None

    minimizer = min((((b[i] / A[i][pivot], i) for i in range(len(b)) if A[i].get(pivot, 0) > 0)), default=None)

    if minimizer == None:
        raise Exception('Unbounded.')

    (_, i0) = minimizer

    return (i0, pivot)

def greedyStrategy(c, A, b, c0, base):
    """
    Greedy naive strategy for pivot choice.

    Returns the pivot (i0,j0). Returns None if the current solution is optimal.
    """
    if next((cval for cval in c.values() if cval > 0), None) == None:
        return None
    
    cchoice = ((c[var], var) for var in c if c[var] > 0)

    minimizers = (min(((b[i] / A[i][var], i, var) for i in range(len(b)) if A[i].get(var, 0) > 0), default=None) for (_,var) in cchoice)

    filtered = (minim for minim in minimizers if minim != None)

    result = max(((v * c[var], i, var) for (v,i,var) in filtered), default=None)

    if result == None:
        raise Exception('Unbounded.')

    (_, i0, pivot) = result
    return (i0, pivot)

###############################
#  Simplex related functions  #
###############################

def newNonConflictingName(baseName, usedNames, *, appendStr='_'):
    """
    Returns a string that is not in usedNames by consecutively appending appendStr to baseName
    until the first is found (zero appends is possible if baseName is not in usedName).
    """
    newName = baseName
    while newName in usedNames:
        newName += appendStr
    return newName

def pivotOperation(c, A, b, c0, base, i0, pivot):
    """
    Applies pivot operation to (c, A, b, c0) removing i0 from base and putting pivot in its place.
    Returns new c0.
    """
    if base[i0] == pivot:
        return c0

    norm = A[i0][pivot]
    varclass = type(norm)

    for var in A[i0]:
        A[i0][var] /= norm

    b[i0] /= norm

    pivotlinevars = set(A[i0])

    for i in range(len(A)):
        if i != i0:
            mult = A[i].get(pivot, 0)
            if mult != 0:
                b[i] -= mult * b[i0]
                for var in pivotlinevars.union(A[i]):
                    Ai0var = A[i0].get(var, 0)
                    if Ai0var != 0:
                        A[i][var] = A[i].get(var, varclass(0)) - mult * Ai0var
                        if A[i][var] == varclass(0):
                            del A[i][var]

    if c == None:
        return None

    objmult = c.get(pivot, 0)
    if objmult != 0:
        for var in pivotlinevars.union(c):
            Ai0var = A[i0].get(var, 0)
            if Ai0var != 0:
                c[var] = c.get(var, varclass(0)) - objmult * Ai0var
                if c[var] == varclass(0):
                    del c[var]
        c0 += objmult * b[i0]
    
    base[i0] = pivot

    return c0

def innerSolveSimplex(c, A, b, *, c0, base, pivotchoice, cutpoint, callback):
    """
    Solves the following linear program
    maximize c0 + c^T x
    s.t. Ax = b
         x >= 0
    by turning it into the canonical tableau of the optimal solution.
    Uses the function pivotchoice(c, A, b, base) to choose the pivot (i,j).

    Returns (x,v) where x is an optimal solution and v is its value
    Raises Exception('Unbounded.') if the problem is unbounded.
    Assumption: x is feasible and base is a list of basic variables;
                the minor of A indexed by basis is an identity;
                c restricted to indices in basis is identically 0

    If cutpoint is given, stops at the first solution with value at least cutpoint.

    Calls callback(**kwargs) on every iteration with the following keys and values.
    - it: iteration number
    - c: objective vector
    - A: system matrix
    - b: system right-hand side
    - c0: current objective value
    - base: base variable indices
    - solved: True if no further iterations are going to be made (it may be due to cutpoint)
    - optimal: True if current point is optimal
    """
    it = 0
    ret = pivotchoice(c, A, b, c0, base)
    while ret and (cutpoint==None or c0 < cutpoint):
        callback(it=it, c=c, A=A, b=b, c0=c0, base=base, solved=False, optimal=False)
        it += 1
        (i0, pivot) = ret
        c0 = pivotOperation(c, A, b, c0, base, i0, pivot)
        ret = pivotchoice(c, A, b, c0, base)

    callback(it=it, c=c, A=A, b=b, c0=c0, base=base, solved=True, optimal=(ret==None))
    x = {}
    for i in range(len(base)):
        x[base[i]] = b[i]
    return (x, c0)

def dot_product(a, b):
    assert len(a) == len(b), 'Length of vectors differ:\n' + str(a) + '\n' + str(b)
    return sum((a[i] * b[i] for i in range(len(a))))

def gramSchmidtFindRankIndices(A, order, *, varnames=None, varclass=None):
    if varnames == None:
        varnames = set()
        for line in A:
            varnames.update(line)

    if varclass == None:
        varclass = type(next(iter(A[0].values())))
    
    columns = dict(((var,[A[i].get(var, varclass(0)) for i in range(len(A))]) for var in varnames))
    base = []
    squarednorms = []

    remainingrank = min(len(A), len(varnames))
    for i in range(len(order)):
        for j in range(len(base)):
            multiplier = dot_product(columns[order[i]], columns[base[j]]) / squarednorms[j]
            for k in range(len(columns[order[i]])):
                columns[order[i]][k] -= multiplier * columns[order[j]][k]

        sn = dot_product(columns[order[i]], columns[order[i]])
        if sn != varclass(0):
            base.append(order[i])
            squarednorms.append(sn)
            remainingrank -= 1
            if remainingrank == 0:
                break
        
    return base

def injectPreProcess(base, candidateBase):
    order = []
    candidateBaseSet = set()
    for v in candidateBase:
        if v not in candidateBaseSet:
            order.append(v)
            candidateBaseSet.add(v)

    return (order, candidateBaseSet)

def injectCandidateBaseViaGS(A, b, base, candidateBase, *, varnames=None, varclass=None, callback):
    """
    Tries to inject candidateBase as a base using Gram--Schmidt without caring about feasibility.

    Calls callback(**kwargs) on every iteration with the following keys and values.
    - it: iteration number
    - c: None
    - A: system matrix
    - b: system right-hand side
    - c0: None
    - base: base variable indices
    """
    (order, candidateBaseSet) = injectPreProcess(base, candidateBase)
    order += [v for v in base if v not in candidateBaseSet]

    GSbaseSet = gramSchmidtFindRankIndices(A, order, varnames=varnames, varclass=varclass)
    remainingBase = [True] * len(base)

    it = 0
    for pivot in GSbaseSet:
        callback(it=it, c=None, A=A, b=b, c0=None, base=base)
        it += 1
        i0 = next((i for i in range(len(A)) if remainingBase[i] and A[i].get(pivot, 0) != varclass(0)), None)
        assert i0 != None, 'injectCandidateBaseViaGS found invalid base.'
        pivotOperation(None, A, b, None, base, i0, pivot)
        remainingBase[i0] = False

def sequentialInjectCandidateBase(c, A, b, c0, base, cutpoint, candidateBase, times, *, callback):
    """
    Tries to inject candidateBase as a base sequentially and keeping feasibility.

    Calls callback(**kwargs) on every iteration with the following keys and values.
    - it: iteration number
    - rep: repetition number
    - c: None
    - A: system matrix
    - b: system right-hand side
    - c0: None
    - base: base variable indices
    """
    (order, _) = injectPreProcess(base, candidateBase)

    t = 0
    while t < times:
        t += 1
        changes = False
        remainingBase = [True] * len(base)

        it = 0
        for pivot in order:
            callback(it=it, rep=t, c=c, A=A, b=b, c0=c0, base=base)
            it += 1

            minvalue = None
            i0 = None
            for i in range(len(A)):
                if A[i].get(pivot, 0) > 0:
                    v = b[i] / A[i][pivot]
                    if minvalue == None or v < minvalue:
                        minvalue = v
                        i0 = i if remainingBase[i] else None
                    elif i0 == None and remainingBase[i] and v == minvalue:
                        i0 = i

            if minvalue == None:
                if c.get(pivot, 0) > 0:
                    raise Exception('Unbounded.')
                continue

            if i0 == None:
                continue
            
            if base[i0] != pivot:
                c0 = pivotOperation(c, A, b, c0, base, i0, pivot)
                if cutpoint != None and c0 >= cutpoint:
                    return c0
                changes = True
            remainingBase[i0] = False

        if changes == False:
            break
    return c0



def solveNamedStandardSimplex(c, A, b, c0=None, *, pivotchoice=basicStrategy, cutpoint=None,
                              candidateBase=None, gramSchmidtInject=False, feasibilityInject=0, solvingInject=0,
                              callback=dummy, createCopy=True):
    """
    Solves the following linear program
    maximize c0 + c^T x
    s.t. Ax <= b
          x >= 0

    Returns a pair (x,v) where x is an optimal solution in the form of a dictionary (absent keys have value 0)
    and v is its value.
    Raises Exception('Infeasible.') if the problem is not feasible.
    Raises Exception('Unbounded.') if the problem is unbounded.

    If candidateBase is given, tries to put its elements as base before solving the problem at the following
    places:
    - If gramSchmidtInject is True, tries to inject before creating auxiliary feasibility problem using
    Gram--Schmidt injection.
    - Tries to inject on the auxiliary feasiblity problem using sequential injection a maximum of
    feasibilityInject times.
    - Tries to inject on the original problem (after feasibility is reached) using sequential
    injection a maximum solvingInject times.

    If cutpoint is given, stops at the first solution with value at least cutpoint.

    Calls callback(**kwargs) on every iteration adding the following keys and values.
    - phase: phase name
    -- callback is also passed to the following functions which may include keys:
    - injectCandidateBaseViaGS
    - sequentialInjectCandidateBase
    - innerSolveSimplex

    If createCopy is True (default), creates a copy of c, A, b before solving. If createCopy is False, c, A, b will
    be destroyed.
    """

    ########## Initial setup ##########
    
    def phaseCallback(phasename):
        return lambda **kwargs : callback(phase=phasename, **kwargs)

    if len(b) > 0:
        varclass = type(b[0])
    elif len(c) > 0:
        elem = next(iter(c.values()), None)
        if elem == None:
            return ({}, 0) # System has no restrictions and has zero objective
        varclass = type(elem)

    if c0 == None:
        c0 = varclass(0)
    if cutpoint != None:
        cutpoint = varclass(cutpoint)

    Ax = A.copy() if createCopy else A
    bx = b.copy() if createCopy else b

    varnames = set(c)

    for line in A:
        varnames.update(line)

    ########## Adding slack variables ##########

    allnames = varnames.copy()
    slacknames = []
    for i in range(len(bx)):
        newname =  newNonConflictingName('_slack%d' % i, allnames)
        slacknames.append(newname)
        allnames.add(newname)
        Ax[i][newname] = varclass(1)

    base = slacknames.copy()

    if candidateBase != None and gramSchmidtInject == True:
        injectCandidateBaseViaGS(Ax, bx, base, candidateBase, callback=phaseCallback('GSInject'),
                                 varnames=allnames, varclass=varclass)

    ########## Feasibility phase ##########
    
    feasibilityName = newNonConflictingName('_feasibility', allnames)
    cfeas = {feasibilityName:varclass(-1)}

    for i in range(len(bx)):
        Ax[i][feasibilityName] = varclass(-1)

    minimizer = min(((b[i],i) for i in range(len(b))), default=None)

    if minimizer != None:
        (value, i0) = minimizer
        if value < varclass(0):
            c0feas = pivotOperation(cfeas, Ax, bx, varclass(0), base, i0, feasibilityName)
            if candidateBase != None and feasibilityInject > 0:
                c0feas = sequentialInjectCandidateBase(c=cfeas, A=Ax, b=bx, c0=c0feas, base=base, cutpoint=None,
                                                       candidateBase=candidateBase, times=feasibilityInject,
                                                       callback=phaseCallback('FeasibilitySeqInject'))
            (x,v) = innerSolveSimplex(cfeas, Ax, bx, c0=c0feas, base=base,
                                      pivotchoice=pivotchoice, cutpoint=None,
                                      callback=phaseCallback('Feasibility'))
            if x.get(feasibilityName, 0) > varclass(0):
                raise Exception('Infeasible.')

    try:
        i0 = base.index(feasibilityName)
        pivot = None
        for (name,val) in Ax[i0].items():
            if val != varclass(0):
                pivot = name
                break
        c0feas = pivotOperation(cfeas, Ax, bx, c0feas, base, i0, pivot)
    except ValueError:
        pass
            
    for i in range(len(bx)):
        Ax[i].pop(feasibilityName, None)

    cx = c.copy() if createCopy else c 

    cbasic = [cx.get(var, varclass(0)) for var in base]

    it = 0
    initialPivotCallback = phaseCallback('InitialPivot')

    for i in range(len(bx)):
        if cbasic[i] != varclass(0):
            initialPivotCallback(it=it, c=cx, A=Ax, b=bx, base=base)
            for var in Ax[i]:
                cx[var] = cx.get(var, varclass(0)) - cbasic[i] * Ax[i][var]
            if cx[var] == varclass(0):
                del cx[var]
            c0 += cbasic[i] * bx[i]
            it += 1

    ########## Solution phase ##########

    if candidateBase != None and solvingInject > 0:
        c0 = sequentialInjectCandidateBase(c=cx, A=Ax, b=bx, c0=c0, base=base, cutpoint=cutpoint,
                                           candidateBase=candidateBase, times=solvingInject,
                                           callback=phaseCallback('SolutionSeqInject'))

    (x,v) = innerSolveSimplex(cx, Ax, bx, c0=c0, base=base, pivotchoice=pivotchoice, 
                              cutpoint=cutpoint,
                              callback=phaseCallback('Solution'))

    for slack in slacknames:
        x.pop(slack, None)
    
    return (x, v)


def dictifySystem(c, A, b, *, removeZeros=False):
    """
    Transforms a system (c, A, b) given in list form to the equivalent system (cx, Ax, bx) given in dict.
    Returns transformed system (cx, Ax, bx).
    If removeZeros is True, removes zero entries from the system.
    """
    if c != None:
        varnames = range(len(c))
    elif A != None and len(A) > 0:
        varnames = range(len(A[0]))

    cx = dict(zip(varnames,c)) if c != None else None
    Ax = [dict(zip(varnames,line)) for line in A] if A != None else None
    bx = b.copy() if b != None else None

    if removeZeros:
        cx = dict((key, val) for (key, val) in cx.items() if val != 0)
        Ax = [dict((key, val) for (key, val) in line.items() if val != 0) for line in Ax]

    return (cx, Ax, bx)

def solveStandardSimplex(c, A, b, c0=None, *, pivotchoice=basicStrategy, cutpoint=None,
                         candidateBase=None, gramSchmidtInject=False, feasibilityInject=0, solvingInject=0,
                         callback=dummy, createCopy=True):
    """
    Wrapper for solveNamedStandardSimplex that transforms a system given in list form to a system given
    in dict form before calling solveNamedStandardSimplex. Always creates a copy of the system beforehand.

    Mostly for backward compatibility.
    """
    (cx, Ax, bx) = dictifySystem(c, A, b)
    (xx, v) = solveNamedStandardSimplex(cx, Ax, bx, c0, pivotchoice=pivotchoice, cutpoint=cutpoint,
                                        candidateBase=candidateBase, gramSchmidtInject=gramSchmidtInject,
                                        feasibilityInject=feasibilityInject, solvingInject=solvingInject,
                                        callback=callback, createCopy=False)
    x = list(range(len(c)))
    for i in range(len(x)):
        x[i] = xx.get(i, 0)

    return (x, v)


def solveNamedLinearSystem(A, b, *, callback=dummy):
    """
    Solves a linear system of the form Ax = b
    by diagonalizing it and returning a vector base of pairs (i,var) such that
    var is the base variable of line i in A.

    The left-hand side is given by A, which is a list of dictionaries.
    The right-hand side is given by b, which is a list of numbers.

    All dictionaries associate variable names to their coefficients.

    Calls callback(**kwargs) on every iteration adding the following keys and values.
    - it: iteration number
    - A: system matrix
    - b: system right-hand side
    - base: (partial) base pairs
    - valid: True if solved system has a solution, False if system does not have a solution, None if system is not
             yet known to have a solution or not.

    Returns a tuple (base, valid), where base is the vector of pairs (i,var) such that
    var is the base variable of line i in A and valid is True if the system has a solution.
    """

    if len(A) != len(b):
        raise Exception('Incompatible system and right-hand side.')

    if len(A) == 0:
        return([], True)

    varclass = type(b[0])

    base = []

    it = 0
    valid = None
    for line in range(len(A)):
        callback(it=it, A=A, b=b, base=base, valid=valid)
        variable = None
        for var in A[line]:
            if A[line][var] != varclass(0):
                variable = var
                break

        if variable == None:
            if b[line] != varclass(0):
                valid = False
            it += 1
            continue

        value = A[line][variable]
        for var in A[line]:
            A[line][var] /= value
        b[line] /= value
        
        for i in range(len(A)):
            if i == line:
                continue
            
            value = A[i].get(variable, None)
            if value == None:
                continue
            for var in A[line]:
                A[i][var] = A[i].get(var, varclass(0)) - value * A[line][var]
            b[i] -= value * b[line]

        base.append((line,variable))
        it += 1

    return (base, valid==None)

def solveNamedSimplex(c, A, b, Aineq, vartype, *, c0=None, maximization=True, varclass=None,
                      pivotchoice=basicStrategy, cutpoint=None, candidateBase=None,
                      gramSchmidtInject=False, feasibilityInject=0, solvingInject=0, callback=dummy):
    """
    Solves a linear program, which is maximization if maximization==True.
    The left-hand side is given by A, which is a list of dictionaries.
    The right-hand side is given by b, which is a list of numbers.
    The objective is given by c, which is a dictionary (plus a constant c0 that defaults to 0).

    All dictionaries associate variable names to their coefficients.
    Aineq[i] gives what type of inequality is given by A[i] and b[i];
        '<=' for A[i] <= b[i], '>=' for A[i] >= b[i] and '=' for A[i] = b[i].
    vartype gives which variables are unbounded ('u'), non-negative ('p') and non-positive ('n').
    varclass (if not None) gives a class of numbers to convert the system to before solving.
    Uses the function pivotchoice(c, A, b, base) to choose the pivot (i,j).

    If candidatebase is given, tries to put its elements as base before solving the problem at the following
    places:
    - If gramSchmidtInject is True, tries to inject before creating auxiliary feasibility problem using
    Gram--Schmidt injection.
    - Tries to inject on the auxiliary feasiblity problem using sequential injection a maximum of
    feasibilityInject times.
    - Tries to inject on the original problem (after feasibility is reached) using sequential
    injection a maximum solvingInject times.

    If cutpoint is given, stops at the first solution with value at least cutpoint if maximization
    and at most cutpoint if minimization.

    Calls callback(**kwargs) on every iteration adding the following keys and values.
    - phase: if in pre-solve phase, this is equal to 'Pre-solve', otherwise this key is not added; it may be 
    added by subfunctions (see below).
    - var2new: dictionary mapping original variable names to new variable names (if variable var split into
    two variables, then var2new[var] will be a pair with both new names, if variable var is pre-solved,
    then var2new[var] will be None).
    - new2var: dictionary mapping new variable names to original variable names
    - maximization: maximization boolean
    - preA: final pre-solve system
    - preb: final pre-solve system right-hand side
    -- callback is also passed to the following functions which may include keys:
    - solveNamedLinearSystem
    - solveNamedStandardSimplex

    Returns a pair (x,v), where v is the optimal value and x is an optimal solution.
    """

    ########## Initial asserts ##########
    
    if len(A) != len(b):
        raise Exception('Incompatible system and right-hand side.')

    if len(A) != len(Aineq):
        raise Exception('Incompatible system and inequality types.')

    ########## Initial setup ##########

    if varclass == None:
        varclass = lambda x : x

    if c0 == None:
        c0 = 0
    c0 = varclass(c0)
    if cutpoint != None:            
        cutpoint = varclass(cutpoint)

    varnames = set(c)
    varnames.update(vartype)

    for line in A:
        varnames.update(line)

    assert candidateBase == None or varnames.issuperset(candidateBase), 'candidateBase not contained in varnames'


    ########## Equality pre-solve ##########
    
    preA = []
    preb = []
    for i in range(len(A)):
        if Aineq[i] == '=':
            preA.append(dict((var,varclass(value)) for (var,value) in  A[i].items()))
            preb.append(varclass(b[i]))

    def preSolveCallback(**kwargs):
        return callback(phase='Pre-solve', **kwargs)

    (presolvebase, valid) = solveNamedLinearSystem(A=preA, b=preb, callback=preSolveCallback)

    if valid == False:
        raise Exception('Infeasible.')

    presolvedvars = dict(((var,i) for (i,var) in presolvebase))

    prec = {}
    for var in varnames:
        prec[var] = varclass(c.get(var,0))
    
    for (i,svar) in presolvebase:
        value = prec[svar]
        if value != 0:
            for v in preA[i]:
                prec[v] -= value * preA[i][v]
            c0 += preb[i] * value

    prec = dict((var, val) for (var, val) in prec.items() if val != 0)

    ########## Variable versus new variable correspondence ##########

    var2new = {}
    new2var = {}
    newnames = set()

    for var in varnames:
        presolveline = presolvedvars.get(var, None)
        if presolveline != None:
            var2new[var] = None
        else:
            if vartype[var] == 'u':
                nameplus = newNonConflictingName(var + '+', newnames)
                newnames.add(nameplus)
                nameminus = newNonConflictingName(var + '-', newnames)
                newnames.add(nameminus)
                var2new[var] = (nameplus, nameminus)
                new2var[nameplus] = new2var[nameminus] = var
            elif vartype[var] == 'p':
                name = newNonConflictingName(var, newnames)
                var2new[var] = name
                new2var[name] = var
                newnames.add(name)
            elif vartype[var] == 'n':
                name = newNonConflictingName(var, newnames)
                var2new[var] = name
                new2var[name] = var
                newnames.add(name)
            else:
                raise Exception("Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var))

    ########## Objective function ##########

    cx = {}

    for var in varnames:
        name = var2new[var]
        if name == None:
            continue

        obj = prec.get(var, None)
        if obj == None:
            continue

        if vartype[var] == 'u':
            (nameplus, nameminus) = name
            cx[nameplus] = obj
            cx[nameminus] = -obj
        elif vartype[var] == 'p':
            cx[name] = obj
        else:
            assert vartype[var] == 'n', "Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var)
            cx[name] = -obj

    ########## New names candidate base ##########

    newNameCandidateBase = None
    if candidateBase != None:
        newNameCandidateBase = []
        for var in candidateBase:
            name = var2new[var]
            if name == None:
                pass
            elif vartype[var] == 'u':
                newNameCandidateBase.extend(name)
            else:
                newNameCandidateBase.append(name)

    ########## Original restrictions ##########

    Ax = []
    bx = []

    for line in range(len(A)):
        if Aineq[line] == '=':
            continue

        rest = {}

        for var in varnames:
            name = var2new[var]
            if name == None:
                continue

            coeff = varclass(A[line].get(var, 0))

            for (i, solvar) in presolvebase:
                Alinesolvar = A[line].get(solvar, 0)
                if Alinesolvar == 0:
                    continue
                Alinesolvar = varclass(Alinesolvar)

                preAivar = preA[i].get(var, 0)
                if preAivar == 0:
                    continue
                preAivar = varclass(preAivar)

                coeff -= Alinesolvar * preAivar

            if coeff == varclass(0):
                continue

            if vartype[var] == 'u':
                (nameplus, nameminus) = name
                rest[nameplus] = coeff
                rest[nameminus] = -coeff
            elif vartype[var] == 'p':
                rest[name] = coeff
            else:
                assert vartype[var] == 'n', "Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var)
                rest[name] = -coeff

        rhs = varclass(b[line])
        for (i, solvar) in presolvebase:
            Alinesolvar = A[line].get(solvar, 0)
            if Alinesolvar == 0:
                continue
            Alinesolvar = varclass(Alinesolvar)

            rhs -= Alinesolvar * preb[i]

        if Aineq[line] == '<=':
            pass
        elif Aineq[line] == '>=':
            rest = dict((name,-val) for (name,val) in rest.items())
            rhs = -rhs
        else:
            raise Exception("Wrong inequality type: ``%s'' for inequality ``%d''." % (Aineq[line],line))

        Ax.append(rest)
        bx.append(rhs)

    ########## Restrictions inherited from pre-solve and solved variable boundedness ##########

    for solvar in presolvedvars:
        line = presolvedvars[solvar]

        if vartype[solvar] == 'u':
            continue

        rest = {}

        for name in newnames:
            var = new2var[name]
            coeff = preA[line].get(var, None)
            if coeff == None:
                continue

            if vartype[var] == 'u':
                (nameplus, nameminus) = var2new[var]
                rest[nameplus] = coeff
                rest[nameminus] = -coeff
            elif vartype[var] == 'p':
                rest[name] = coeff
            else:
                assert vartype[var] == 'n', "Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var)
                rest[name] = -coeff

        rhs = preb[line]
        
        if vartype[solvar] == 'p':
            pass
        elif vartype[solvar] == 'n':
            rest = dict((name,-val) for (name,val) in rest.items())
            rhs = -rhs
        else:
            raise Exception("Wrong variable type: ``%s'' for ``%s''." % (vartype[solvar], solvar))

        Ax.append(rest)
        bx.append(rhs)

    ########## Maximization/minimization ##########

    if maximization == False:
        cx = dict((name,-val) for (name,val) in cx.items())
        c0 = -c0
        if cutpoint != None:
            cutpoint = -cutpoint
    
    ########## Solver function call ##########

    def namedCallback(**kwargs):
        return callback(var2new=var2new, new2var=new2var, maximization=maximization,
                        preA=preA, preb=preb, **kwargs)

    (xx, value) = solveNamedStandardSimplex(c=cx, A=Ax, b=bx, c0=c0, pivotchoice=pivotchoice, cutpoint=cutpoint,
                                            candidateBase=newNameCandidateBase,
                                            gramSchmidtInject=gramSchmidtInject,
                                            feasibilityInject=feasibilityInject, solvingInject=solvingInject,
                                            callback=namedCallback, createCopy=False)

    ########## Solution translation ##########
                    
    if maximization == False:
        value = -value

    x = {}

    for var in varnames:
        presolveline = presolvedvars.get(var, None)
        if presolveline == None:
            if vartype[var] == 'u':
                x[var] = xx.get(var2new[var][0], varclass(0)) - xx.get(var2new[var][1], varclass(0))
            elif vartype[var] == 'p':
                x[var] = xx.get(var2new[var], varclass(0))
            else:
                assert vartype[var] == 'n', "Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var)
                x[var] = -xx.get(var2new[var], varclass(0))

    for solvar in presolvedvars:
        line = presolvedvars[solvar]
        x[solvar] = preb[line]
        for var in preA[line]:
            if var != solvar:
                x[solvar] -= preA[line][var] * x.get(var, varclass(0))
                
    return (x, value)

def checkFeasibility(A, b, Aineq, vartype, x):
    """
    Checks whether x is a feasible point of a linear program.
    The left-hand side is given by A, which is a list of dictionaries.
    The right-hand side is given by b, which is a list of numbers.

    All dictionaries associate variable names to their coefficients.
    Aineq[i] gives what type of inequality is given by A[i] and b[i];
        '<=' for A[i] <= b[i], '>=' for A[i] >= b[i] and '=' for A[i] = b[i].
    vartype gives which variables are unbounded ('u'), non-negative ('p') and non-positive ('n').

    Returns a pair (feasible, slack), where feasible is True if and only if x is feasible and slack is a list of
    slacknesses of the system that are defined as the left-hand side minus the right-hand side, that is, we have
    slack[i] = <A[i], x> - b[i].
    
    The solution is feasible if and only if x respects vartype and the following hold for every i.
    - If Aineq[i] = '<=', then slack[i] <= 0.
    - If Aineq[i] = '=', then slack[i] = 0.
    - If Aineq[i] = '>=', then slack[i] >= 0.
    """

    ########## Initial asserts ##########
    
    if len(A) != len(b):
        raise Exception('Incompatible system and right-hand side.')

    if len(A) != len(Aineq):
        raise Exception('Incompatible system and inequality types.')

    ########## Slack computation ##########

    slack = [sum(A[i].get(varname, 0) * x[varname] for varname in x) - b[i] for i in range(len(b))]

    ########## Feasibility check ##########

    feasible = True

    for var in x:
        if vartype[var] == 'u':
            pass
        elif vartype[var] == 'p':
            if x[var] < 0:
                feasible = False
        elif vartype[var] == 'n':
            if x[var] > 0:
                feasible = False
        else:
            raise Exception("Wrong variable type: ``%s'' for ``%s''." % (vartype[var], var))

    for i in range(len(A)):
        if Aineq[i] == '=':
            if slack[i] != 0:
                feasible = False
        elif Aineq[i] == '<=':
            if slack[i] > 0:
                feasible = False
        elif Aineq[i] == '>=':
            if slack[i] < 0:
                feasible = False
        else:
            raise Exception("Wrong inequality type: ``%s'' for inequality ``%d''." % (Aineq[i],i))

    return (feasible, slack)

def checkValue(c, x, v, *, c0=0):
    """
    Checks whether x gives objective value v for a given linear program whose objective is c.

    Returns a pair (correct, truevalue), where truevalue = <c, x> + c0 and correct is v==truevalue.
    """

    truevalue = sum(c.get(var, 0) * x[var] for var in x) + c0

    return (v == truevalue, truevalue)
