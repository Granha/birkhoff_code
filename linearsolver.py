import Partition
import simplex
import fractions
import pickle
import math
import glob
import re
import os
import multiprocessing
import gc
import CompactSet
import decimal

#######################################
#  Rounding functions and generators  #
#######################################

def dummy(rest, rhs, **kwargs):
    """
    Returns rhs
    """
    return rhs

def joinPrimalRounders(*args):
    """
    Receives a list of primal rounders and returns a primal rounder that executes them in order.
    """
    def joined(rest, rhs, **kwargs):
        for f in args:
            rhs = f(rest, rhs, **kwargs)
    return rhs

def generateMildKostkaRounder(cutfraction=fractions.Fraction(1,100), rescale=True):
    """
    Returns a function that performs mild (primal) rounding with parameter cutfraction.

    If rescale is True, then rescales inequalities.
    """
    def rounder(rest, rhs, **kwargs):
        cutpoint = cutfraction * fractions.Fraction(rhs)
        divisor = min((v for v in rest.values() if v >= cutpoint), default=None)
        if divisor == None:
            raise Exception('Cutpoint %s is too large for restriction.\n' % str(simplex.toDecimal(cutpoint)) +
                            'Largest coefficient is %s\n' % str(simplex.toDecimal(max(v for v in rest.values()))))

        divisor = fractions.Fraction(divisor)
        
        for i in rest:
            rest[i] = math.floor(fractions.Fraction(rest[i]) / divisor)
            if rescale==False:
                rest[i] *= divisor

        if rescale:
            return math.ceil(fractions.Fraction(rhs) / divisor)
        else:
            return math.ceil(fractions.Fraction(rhs) / divisor) * divisor

    return rounder

def generateAggressiveKostkaRounder(radius=20):
    """
    Returns a function that performs aggressive (primal) rounding with parameter radius.
    """
    assert radius > 0, 'Radius is negative.'

    def rounder(rest, rhs, **kwargs):
        pairs = [(v,name) for (name,v) in rest.items() if v > 0]
        pairs.sort()

        for (v,name) in pairs[:-radius]:
            rest[name] = type(rest[name])(0)

        divisor,_ = pairs[-radius] if radius <= len(pairs) else pairs[0]
        divisor = fractions.Fraction(divisor)

        for i in rest:
            rest[i] = math.floor(fractions.Fraction(rest[i]) / divisor)
        return fractions.Fraction(rhs) / divisor

    return rounder

def dualDummy(obj, A, b, Aineq, vartype, obj0, **kwargs):
    """
    Returns obj0
    """
    return obj0

def joinDualRounders(*args):
    """
    Receives a list of dual rounders and returns a dual rounder that executes them in order.
    """
    def joined(obj, A, b, Aineq, vartype, obj0, **kwargs):
        for f in args:
            obj0 = f(obj, A, b, Aineq, vartype, obj0, **kwargs)
    return obj0


def generateMildDualKostkaRounder(cutpoint=fractions.Fraction(1,100)):
    """
    Returns a function that performs mild dual rounding with parameter cutpoint.
    """
    def rounder(obj, A, b, Aineq, vartype, obj0, **kwargs):
        for i in range(len(A)):
            if Aineq[i] == '>=':
                for varname in A[i]:
                    if varname[0] == 'y':
                        A[i][varname] /= abs(obj[varname])
                        if A[i][varname] < cutpoint:
                            A[i][varname] = type(A[i][varname])(0)
                        
        for varname in obj:
            if varname[0] == 'y':
                obj[varname] /= abs(obj[varname])

        return obj0

    return rounder

def intlog(val, base, *, floor=True):
    """
    Returns a pair (lg, n), where lg is the either the floor (if floor=True) or the ceiling (if floor=False)
    of the logarithm of val on base base and n = base ** lg.

    Assumes that val > 0 and base > 0 and base != 1.
    """
    sgn = 1
    flipval = False
    if base < 1:
        base = 1/base
        sgn = -sgn
    if val < 1:
        val = 1/val
        flipval = True
        sgn = -sgn
    lg = 0
    lastn = 1
    n = 1
    while val > n:
        lastn = n
        n = lastn * base
        lg += 1
    if val < n and ((sgn == -1) != floor):
        lg -= 1
        n = lastn
    if flipval:
        n = 1/n
    return (sgn * lg, n)

def generateAggressiveDualKostkaRounder(cutlog=None, roundObjective=True):
    """
    Returns a function that performs aggressive dual rounding with parameter cutpoint.
    """
    def rounder(obj, A, b, Aineq, vartype, obj0, **kwargs):
        c = kwargs['c']

        if roundObjective:
            for varname in obj:
                if varname[0] == 'y':
                    (lg, n) = intlog(obj[varname], c)
                    obj[varname] = n

        for i in range(len(A)):
            if Aineq[i] == '>=':
                for varname in A[i]:
                    if varname[0] == 'y':
                        (lg, n) = intlog(A[i][varname], c)
                        if cutlog != None and lg < cutlog:
                            A[i][varname] = type(A[i][varname])(0)
                        else:
                            A[i][varname] = n / abs(obj[varname])
                        
        for varname in obj:
            if varname[0] == 'y':
                obj[varname] /= abs(obj[varname])

        return obj0

    return rounder

###################
#  Parametrizers  #
###################

def dummyParametrizer():
    """
    Dummy function used for default value.
    """
    return

def generateExponentialParametrizer(base, initial=0):
    """
    Returns a parametrizer (for either Y or W, see createBirkhoffDual) f(m) that returns:
    - None, if m < initial.
    - base ** m, if m >= inital.
    """
    return (lambda m: None if m < initial else base ** m)


###################################################
#  Birkhoff problem creators, solvers and savers  #
###################################################

def tailbound(l, k0, c):
    """
    Returns the tailbound associated with Parseval level l with truncation parameter k0 and
    pseudorandomness parameter c.
    """
    return fractions.Fraction(15,10) * ((2*(k0 + l) + 4) * (c/2)**(2*(k0 + l) + 4) - (2*(k0 + l)) * (c/2)**(2*(k0 + l) + 8))/(1 - (c/2)**4)**2

def maxRestrictionsBirkhoffPrimal(lmax, *, lset=None):
    """
    Returns a generator that generates the (non-rounded) max restrictions of the primal linear program
    associated to the Birkhoff problem.

    Each item returned by the generator is a tuple (lhs, ineq, rhs), where
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    The parameters of this function are:
    - lmax is the maximum Parseval level used.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    """
    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    for l in lset:
        psiname = 'psi%03d' % l
        yield ({psiname:1, 'Max':-1}, '<=', 0)

def parsevalRestrictionsBirkhoffPrimal(lmax, k0 = 15, c = fractions.Fraction(17,10), *, lset=None,
                                       negativeBelly=False, negativeLeg=False,
                                       precomputedPartitionList=None):
    """
    Returns a generator that generates the (non-rounded) Parseval restrictions of the primal linear program
    associated to the Birkhoff problem.

    Each item returned by the generator is a tuple (lhs, ineq, rhs), where
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    The parameters of this function are:
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.
    """
    partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    for l in lset:
        parsevaltable = Partition.parseval(l, k0, negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                           precomputedPartitionList=partitions[l])
        psiname = 'psi%03d' % l

        rest = {}
        rest[psiname] = -1

        for bellysize in range(0, l+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    key = part.bname()
                    rest['%s[%03d]' % (key, k)] = 2 * parsevaltable[k].get(key, 0)

        yield (rest, '=', -2 + 2 * tailbound(l, k0, c))


def youngRestrictionsBirkhoffPrimal(lmax, k0 = 15, c = fractions.Fraction(17,10), *, mset=None,
                                    precomputedPartitionList=None):
    """
    Returns a generator that generates the (non-rounded) Young restrictions of the primal linear program
    associated to the Birkhoff problem.

    Each item returned by the generator is a tuple (lhs, ineq, rhs), where
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    The parameters of this function are:
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.
    """
    partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList
    mmax = 2*k0 + 2*lmax

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    for m in mset:
        rest = {}

        for bellysize in range(0,lmax+1):
            for belly in partitions[bellysize]:
                for k in range(max(belly.height(),1), k0+1):
                    key = belly.bname()
                    rest['%s[%03d]' % (key, k)] = Partition.kostkaNumber(m, k, belly)

        yield (rest, '<=', c ** m - 1)


def createBirkhoffPrimal(lmax, k0 = 15, c = fractions.Fraction(17,10), *, lset=None, mset=None,
                         kostkaRounder=dummy, negativeBelly=False, negativeLeg=False,
                         precomputedPartitionList=None):
    """
    Creates the primal linear program associated to the Birkhoff problem.
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - kostkaRounder is a function to perform primal rounding of Kostka coefficients.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.

    Returns a tuple (obj, A, b, Aineq, vartype, obj0, maximization) to be used as parameters for
    simplex.solveNamedSimplex respectively as (c, A, b, Aineq, vartype, c0, maximization).
    """
    partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList
    mmax = 2*k0 + 2*lmax

    A = []
    b = []
    obj = {}
    Aineq = []
    vartype = {}

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    obj['Max'] = 1
    vartype['Max'] = 'u'

    for l in lset:
        psiname = 'psi%03d' % l
        vartype[psiname] = 'u'

    for bellysize in range(0,lmax+1):
        for belly in partitions[bellysize]:
            for k in range(max(belly.height(),1),k0+1):
                vartype['%s[%03d]' % (belly.bname(), k)] = 'p'

    # Max restrictions
    for (lhs, ineq, rhs) in maxRestrictionsBirkhoffPrimal(lmax=lmax, lset=lset):
        A.append(lhs)
        Aineq.append(ineq)
        b.append(rhs)

    # Parseval restrictions
    for (lhs, ineq, rhs) in parsevalRestrictionsBirkhoffPrimal(lmax=lmax, k0=k0, c=c, lset=lset,
                                                               negativeBelly=negativeBelly,
                                                               negativeLeg=negativeLeg,
                                                               precomputedPartitionList=partitions):
        A.append(lhs)
        Aineq.append(ineq)
        b.append(rhs)

    # Young restrictions
    for (lhs, ineq, rhs) in youngRestrictionsBirkhoffPrimal(lmax=lmax, k0=k0, c=c, mset=mset,
                                                            precomputedPartitionList=partitions):
        # Primal rounding on-the-fly
        rhs = kostkaRounder(lhs, rhs)

        A.append(lhs)
        Aineq.append(ineq)
        b.append(rhs)

    return (obj, A, b, Aineq, vartype, 0, False)

def solveBirkhoffPrimal(lmax, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                        pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                        lset=None, mset=None, kostkaRounder=dummy, negativeBelly=False, negativeLeg=False,
                        precomputedPartitionList=None):
    """
    Creates and solves the primal linear program associated to the Birkhoff problem using the functions
    createBirkhoffPrimal and simplex.solveNamedSimplex.

    See createBirkhoffPrimal for meaning of parameters lmax, k0, c, lset, mset, kostkaRounder,
    negativeBelly, negativeLeg, precomputedPartitionList.

    See simplex.solveNamedSimplex for meaning of parameters varclass, pivotchoice, cutpoint, callback.

    Returns the result of the latter function, that is, 
    returns (x,v), where x is an optimal solution and v is an optimal value.
    """
    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffPrimal(lmax, k0, c,
                                                                       lset=lset, mset=mset,
                                                                       kostkaRounder=kostkaRounder,
                                                                       negativeBelly=negativeBelly,
                                                                       negativeLeg=negativeLeg,
                                                                       precomputedPartitionList=precomputedPartitionList)

    return simplex.solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=c0,
                                     maximization=maximization, varclass=varclass,
                                     pivotchoice=pivotchoice, cutpoint=cutpoint, callback=callback)

def distributionRestrictionsBirkhoffDual(lmax, *, lset=None):
    """
    Returns a generator that generates the (non-rounded) distribution restriction of the dual linear program
    associated to the Birkhoff problem.

    Each item returned by the generator is a tuple (lhs, ineq, rhs), where
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    The parameters of this function are:
    - lmax is the maximum Parseval level used.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    """
    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    yield (dict([('w[%03d]' % l, 1) for l in lset]), '=', 1)

def generateSinglePartitionRestriction(k, belly, lset, mset,
                                       negativeBelly=False, negativeLeg=False):
    """
    Returns the partition restriction of the dual linear program associated to the Birkhoff problem
    corresponding to belly belly and leg k.

    The restriction is in the form (k, belly, lhs, ineq, rhs), where
    - k is the leg of the partition to which the restriction is associated.
    - belly is the belly of the partition to which the restriction is associated.
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    See also partitionRestrictionsBirkhoffDual.

    The parameters of this function are:
    - lset is the set of all Parseval levels used.
    - mset is the set of all Young restrictions used.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    """
    rest = {('w[%03d]' % l): 2 * Partition.singlePartitionParseval(l, k, belly,
                                                                   negativeBelly=negativeBelly,
                                                                   negativeLeg=negativeLeg)
            for l in lset}

    for m in mset:
        rest['y[%03d]' % m] = Partition.kostkaNumber(m, k, belly)

    return (k, belly, rest, '>=', 0)

def generateSingleLargeLegPartitionRestriction(k, belly, lset, mset, negativeLeg=False):
    """
    Returns the partition restriction of the dual linear program associated to the Birkhoff problem
    corresponding to bellies represented by JointPartition belly and leg k.

    The restriction is in the form (k, belly, lhs, ineq, rhs), where
    - k is the leg of the partition to which the restriction is associated.
    - belly is the JointPartition representing bellies of the partitions to which the restriction
    is associated.
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    See also partitionRestrictionsBirkhoffDual.

    The parameters of this function are:
    - lset is the set of all Parseval levels used.
    - mset is the set of all Young restrictions used.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    """
    rest = {'w[%03d]' % belly.size():(-2 if (k % 2 == 1 or negativeLeg) else 2)}

    for m in mset:
        rest['y[%03d]' % m] = Partition.kostkaBinomial(m, k, belly.size()) * Partition.binomial(k, belly.size())
        # The extra factor comes from
        # f_\mu \geq \binom{k}{\lvert\beta\rvert} f_\beta,
        # where \mu is the partition that appears in the Kostka computation.

    return (k, belly, rest, '>=', 0)

def partitionRestrictionsBirkhoffDual(lmax, k0 = 15, *, lset=None, mset=None,
                                      negativeBelly=False, negativeLeg=False, restrictionList=None,
                                      jointLargeLeg=False, precomputedPartitionList=None):
    """
    Returns a generator that generates the (non-rounded) partition restrictions of the dual linear program
    associated to the Birkhoff problem.

    Each item returned by the generator is a tuple (k, belly, lhs, ineq, rhs), where
    - k is the leg of the partition to which the restriction is associated.
    - belly is the belly of the partition to which the restriction is associated.
    - lhs is a dictionary associating variables to coefficients of the left-hand side;
    - ineq is either '<=', '>=' or '=' indicating the type of inequality;
    - rhs is the right-hand side.

    The parameters of this function are:
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - restrictionList is a list (or generator) of which partitions should have their restrictions included. Each
    partition in this list is encoded as a tuple (k, belly), where k is its leg and belly is its belly. If this 
    parameter is None (default), then the function includes all partitions of belly with size at most lmax and leg
    at most k0.
    - if jointLargeLeg is True, replaces all restrictions of corresponding to partitions (k, belly) with large
    leg (i.e., k > lmax) and same belly size s by a single restriction (k, jp, lhs, ineq, rhs)
    that implies them all, where jp is an instance of the class Partition.JointPartition such that jp.size() = s.
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function (this computation only takes place
    if restrictionList is None).
    """
    mmax = 2*k0 + 2*lmax

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    if restrictionList == None:
        partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList
        restrictionList = ((k, belly) for bellysize in range(0, lmax+1) for belly in partitions[bellysize] for k in range(max(belly.height(), 1), k0+1))

    largeLegPairs = set()
        
    for (k, belly) in restrictionList:
        if k > lmax and jointLargeLeg:
            largeLegPairs.add((k, Partition.JointPartition(belly.size())))
            continue

        yield generateSinglePartitionRestriction(k=k, belly=belly, lset=lset, mset=mset,
                                                 negativeBelly=negativeBelly,
                                                 negativeLeg=negativeLeg)

    for (k, jp) in largeLegPairs: # if jointLargeLeg == False, this is empty
        if jp.size() in lset and (k % 2 == 1 or negativeLeg):
            # Adds only if Parseval would give some negative coefficient.
            # (If all coefficients are non-negative, then the restriction is trivially satisfied.)
            
            yield generateSingleLargeLegPartitionRestriction(k=k, belly=jp, lset=lset, mset=mset,
                                                             negativeLeg=negativeLeg)


def createBirkhoffDual(lmax, k0 = 15, c = fractions.Fraction(17,10), *, lset=None, mset=None,
                       kostkaRounder=dummy, dualKostkaRounder=dualDummy, negativeBelly=False, negativeLeg=False,
                       restrictionList=None, jointLargeLeg=False, precomputedPartitionList=None,
                       Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer):
    """
    Creates the dual linear program associated to the Birkhoff problem.
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - kostkaRounder is a function to perform (dualization of the) primal rounding of Kostka coefficients.
    - dualKostkaRounder is a function to perform dual rounding of Kostka coefficients.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - restrictionList is a list (or generator) of which partitions should have their restrictions included. Each
    partition in this list is encoded as a tuple (k, belly), where k is its leg and belly is its belly. If this 
    parameter is None (default), then the function includes all partitions of belly with size at most lmax and leg
    at most k0.
    - if jointLargeLeg is True, replaces all restrictions of corresponding to partitions (k, belly) with large
    leg (i.e., k > lmax) and same belly size s by a single restriction (k, s, lhs, ineq, rhs)
    that implies them all.
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.
    - Yparametrizer: if not dummyParametrizer, must be a fuction that parametrizes each y[m] variable according
    to the returned value of Yparametrizer(m). If Yparametrizer(m) == None, then y[m] is a free variable;
    else y[m] is made to satisfy
    y[m] = Yparametrizer(m) * Y + Yc,
    where Y and Yc are new variables (the same for all m).
    - Wparametrizer: if not dummyParametrizer, must be a fuction that parametrizes each w[l] variable according
    to the returned value of Wparametrizer(l). If Wparametrizer(l) == None, then w[l] is a free variable;
    else w[l] is made to satisfy
    w[l] = Wparametrizer(l) * W + Wc,
    where W and Wc are a new variables (the same for all l).

    WARNING: it is not a priori clear what is the theoretical meaning of applying (dualization of) primal 
    rounding or dual rounding when the parameter jointLargeLeg is True.

    Returns a tuple (obj, A, b, Aineq, vartype, obj0, maximization) to be used as parameters for
    simplex.solveNamedSimplex respectively as (c, A, b, Aineq, vartype, c0, maximization).
    """

    mmax = 2*k0 + 2*lmax

    A = []
    b = []
    obj = {}
    Aineq = []
    vartype = {}

    obj0 = 2

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    # Objective value
    for l in lset:
        obj['w[%03d]' % l] = -2 * tailbound(l, k0, c)

        vartype['w[%03d]' % l] = 'p'

    for m in mset:
        vartype['y[%03d]' % m] = 'p'
        obj['y[%03d]' % m] = - (c ** m - 1) # This may be changed later by primal rounding

    # Distribution restriction
    for (lhs, ineq, rhs) in distributionRestrictionsBirkhoffDual(lmax=lmax, lset=lset):
        A.append(lhs)
        Aineq.append(ineq)
        b.append(rhs)

    # Partition restrictions are computed here but not added to allow for primal rounding
    partRests = partitionRestrictionsBirkhoffDual(lmax=lmax, k0=k0, lset=lset, mset=mset,
                                                  negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                                  restrictionList=restrictionList,
                                                  jointLargeLeg=jointLargeLeg,
                                                  precomputedPartitionList=precomputedPartitionList)

    # Primal rounding
    if kostkaRounder != dummy:
        partRests = list(partRests)

        for m in mset:
            primalrest = {}
            for (k, belly, lhs, ineq, rhs) in partRests:
                primalrest[(k, belly.bname())] = lhs['y[%03d]' % m]
            primalrhs = kostkaRounder(primalrest, -obj['y[%03d]' % m])

            obj['y[%03d]' % m] = -primalrhs
            for (k, belly, lhs, ineq, rhs) in partRests:
                lhs['y[%03d]' % m] = primalrest[(k, belly.bname())]

    # Adding partition restrictions to system
    for (k, belly, lhs, ineq, rhs) in partRests:
        A.append(lhs)
        Aineq.append(ineq)
        b.append(rhs)

    # Dual rounding
    obj0 = dualKostkaRounder(obj, A, b, Aineq, vartype, obj0, c=c)

    # Y parametrization
    if Yparametrizer != dummyParametrizer:
        vartype['Y'] = 'p'
        vartype['Yc'] = 'u'
        for m in mset:
            Yp = Yparametrizer(m)
            if Yp != None:
                A.append({('y[%03d]' % m): -1, 'Y':Yp, 'Yc':1})
                Aineq.append('=')
                b.append(0)

    # W parametrization
    if Wparametrizer != dummyParametrizer:
        vartype['W'] = 'p'
        vartype['Wc'] = 'u'
        for l in lset:
            Wp = Wparametrizer(l)
            if Wp != None:
                A.append({('w[%03d]' % l): -1, 'W':Wp, 'Wc':1})
                Aineq.append('=')
                b.append(0)

    return (obj, A, b, Aineq, vartype, obj0, True)


def solveBirkhoffDual(lmax, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                      pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                      lset=None, mset=None, kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                      negativeBelly=False, negativeLeg=False, restrictionList=None,
                      jointLargeLeg=False, precomputedPartitionList=None,
                      Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer):
    """
    Creates and solves the dual linear program associated to the Birkhoff problem using the functions
    createBirkhoffDual and simplex.solveNamedSimplex.

    See createBirkhoffDual for meaning of parameters lmax, k0, c, lset, mset, kostkaRounder, dualKostkaRounder,
    negativeBelly, negativeLeg, jointLargeLeg, precomputedPartitionList, Yparametrizer, Wparametrizer.

    See simplex.solveNamedSimplex for meaning of parameters varclass, pivotchoice, cutpoint, callback.

    Returns the result of the latter function, that is, 
    returns (x,v), where x is an optimal solution and v is an optimal value.
    """

    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffDual(lmax, k0, c, lset=lset, mset=mset,
                                                                     kostkaRounder=kostkaRounder,
                                                                     dualKostkaRounder=dualKostkaRounder,
                                                                     negativeBelly=negativeBelly,
                                                                     negativeLeg=negativeLeg,
                                                                     restrictionList=restrictionList,
                                                                     jointLargeLeg=jointLargeLeg,
                                                                     precomputedPartitionList=precomputedPartitionList,
                                                                     Yparametrizer=Yparametrizer,
                                                                     Wparametrizer=Wparametrizer)
    return simplex.solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=c0,
                                     maximization=maximization, varclass=varclass,
                                     pivotchoice=pivotchoice, cutpoint=cutpoint, callback=callback)

def solveBirkhoffBoth(lmax, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                      pivotchoice=simplex.greedyStrategy, primalcutpoint=None, dualcutpoint=None,
                      callback=simplex.dummy,
                      lset=None, mset=None, kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                      negativeBelly=False, negativeLeg=False, restrictionList=None,
                      jointLargeLeg=False, precomputedPartitionList=None,
                      Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer):
    """
    Creates and solves both the primal and dual linear programs associated to the Birkhoff problem using
    the functions solveBirkhoffPrimal and solveBirkhoffDual. See these functions for meanings of parameters 
    with the same name.
    - primalcutpoint is cutpoint used with solveBirkhoffPrimal.
    - dualcutpoint is cutpoint used with solveBirkhoffDual.

    Returns a tuple (xp, vp, xd, vd), where (xp, vp) is the return value of the first function and (xd, vd) is
    the return value of the second function.
    """
    (xp, vp) = solveBirkhoffPrimal(lmax=lmax, k0=k0, c=c, varclass=varclass,
                                   pivotchoice=pivotchoice, cutpoint=primalcutpoint, callback=callback,
                                   lset=lset, mset=mset, kostkaRounder=kostkaRounder,
                                   negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                   precomputedPartitionList=precomputedPartitionList)
    (xd, vd) = solveBirkhoffDual(lmax=lmax, k0=k0, c=c, varclass=varclass, 
                                 pivotchoice=pivotchoice, cutpoint=dualcutpoint, callback=callback,
                                 lset=lset, mset=mset,
                                 kostkaRounder=kostkaRounder, dualKostkaRounder=dualKostkaRounder,
                                 negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                 restrictionList=restrictionList, jointLargeLeg=jointLargeLeg,
                                 precomputedPartitionList=precomputedPartitionList,
                                 Yparametrizer=Yparametrizer,
                                 Wparametrizer=Wparametrizer)
    return (xp, vp, xd, vd)


def solveBirkhoffPrimalAndSave(problemfilename=None, solutionfilename=None,
                               lmax=2, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                               pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                               lset=None, mset=None, kostkaRounder=dummy, fileprefix='',
                               negativeBelly=False, negativeLeg=False,
                               precomputedPartitionList=None):
    """
    Creates, solves and saves (both the problem and solution) the primal linear program associated to the Birkhoff
    problem using the functions createBirkhoffPrimal and simplex.solveNamedSimplex.

    See createBirkhoffPrimal for meaning of parameters lmax, k0, c, lset, mset, kostkaRounder,
    negativeBelly, negativeLeg, precomputedPartitionList.

    See simplex.solveNamedSimplex for meaning of parameters varclass, pivotchoice, cutpoint, callback.

    The other parameters are:
    - problemfilename: name of file to save the problem (if None, a standardized name will be generated).
    - solutionfilename: name of file to save the solution (if None, a standardized name will be generated).
    - fileprefix: string to be prepended to standardized names (if any is generated).

    The problem and solution are saved on file as dictionaries using the pickle module.

    Returns the result of the latter function, that is, 
    returns (x,v), where x is an optimal solution and v is an optimal value.
    """
    lsetstr = '' if lset==None else '_lset'
    msetstr = '' if mset==None else '_mset'
    if problemfilename == None:
        problemfilename = fileprefix + 'primalprob_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'
    if solutionfilename == None:
        solutionfilename = fileprefix + 'primalsol_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'

    parameters = {
        'problemfilename':problemfilename,
        'solutionfilename':solutionfilename,
        'lmax':lmax,
        'k0':k0,
        'c':c,
        'varclass':varclass,
        'pivotchoice':pivotchoice.__name__,
        'cutpoint':cutpoint,
        'callback':callback.__name__,
        'lset':lset,
        'mset':mset,
        'kostkaRounder':kostkaRounder.__name__,
        'fileprefix':fileprefix,
        'negativeBelly':negativeBelly,
        'negativeLeg':negativeLeg,
        'precomputedPartitionList':precomputedPartitionList
    }
        
    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffPrimal(lmax, k0, c,
                                                                       lset=lset, mset=mset,
                                                                       kostkaRounder=kostkaRounder,
                                                                       negativeBelly=negativeBelly,
                                                                       negativeLeg=negativeLeg,
                                                                       precomputedPartitionList=precomputedPartitionList)

    
    problem = {'c': c,
               'A': A,
               'b': b,
               'Aineq': Aineq,
               'vartype': vartype,
               'c0': c0,
               'maximization': maximization}

    d = {'type':'primal', 'parameters':parameters, 'problem':problem}

    with open(problemfilename, 'wb') as f:
        pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)

    solution = simplex.solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=c0,
                                         maximization=maximization, varclass=varclass,
                                         pivotchoice=pivotchoice, cutpoint=cutpoint, callback=callback)

    sold = {'type':'primal', 'parameters':parameters, 'solution':solution}

    with open(solutionfilename, 'wb') as f:
        pickle.dump(sold, f, pickle.HIGHEST_PROTOCOL)

    return solution

def solveBirkhoffDualAndSave(problemfilename=None, solutionfilename=None,
                             lmax=2, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                             pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                             lset=None, mset=None, fileprefix='',
                             kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                             negativeBelly=False, negativeLeg=False, restrictionList=None,
                             jointLargeLeg=False, precomputedPartitionList=None,
                             Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer):
    """
    Creates, solves and saves (both the problem and solution) the dual linear program associated to the Birkhoff
    problem using the functions createBirkhoffDual and simplex.solveNamedSimplex.

    See createBirkhoffDual for meaning of parameters lmax, k0, c, lset, mset, kostkaRounder, dualKostkaRounder,
    negativeBelly, negativeLeg, jointLargeLeg, precomputedPartitionList, Yparametrizer, Wparametrizer.

    See simplex.solveNamedSimplex for meaning of parameters varclass, pivotchoice, cutpoint, callback.

    The other parameters are:
    - problemfilename: name of file to save the problem (if None, a standardized name will be generated).
    - solutionfilename: name of file to save the solution (if None, a standardized name will be generated).
    - fileprefix: string to be prepended to standardized names (if any is generated).

    The problem and solution are saved on file as dictionaries using the pickle module.

    Returns the result of the latter function, that is, 
    returns (x,v), where x is an optimal solution and v is an optimal value.
    """
    lsetstr = '' if lset==None else '_lset'
    msetstr = '' if mset==None else '_mset'
    if problemfilename == None:
        problemfilename = fileprefix + 'dualprob_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'
    if solutionfilename == None:
        solutionfilename = fileprefix + 'dualsol_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'

    parameters = {
        'problemfilename':problemfilename,
        'solutionfilename':solutionfilename,
        'lmax':lmax,
        'k0':k0,
        'c':c,
        'varclass':varclass,
        'pivotchoice':pivotchoice.__name__,
        'cutpoint':cutpoint,
        'callback':callback.__name__,
        'lset':lset,
        'mset':mset,
        'fileprefix':fileprefix,
        'kostkaRounder':kostkaRounder.__name__,
        'dualKostkaRounder':dualKostkaRounder.__name__,
        'negativeBelly':negativeBelly,
        'negativeLeg':negativeLeg,
        'restrictionList':restrictionList,
        'jointLargeLeg':jointLargeLeg,
        'precomputedPartitionList':precomputedPartitionList,
        'Yparametrizer':Yparametrizer.__name__,
        'Wparametrizer':Wparametrizer.__name__
    }

    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffDual(lmax, k0, c, lset=lset, mset=mset,
                                                                     kostkaRounder=kostkaRounder,
                                                                     dualKostkaRounder=dualKostkaRounder,
                                                                     negativeBelly=negativeBelly,
                                                                     negativeLeg=negativeLeg,
                                                                     restrictionList=restrictionList,
                                                                     jointLargeLeg=jointLargeLeg,
                                                                     precomputedPartitionList=precomputedPartitionList,
                                                                     Yparametrizer=Yparametrizer,
                                                                     Wparametrizer=Wparametrizer)
    problem = {'c': c,
               'A': A,
               'b': b,
               'Aineq': Aineq,
               'vartype': vartype,
               'c0': c0,
               'maximization': maximization}

    d = {'type':'dual', 'parameters':parameters, 'problem':problem}

    with open(problemfilename, 'wb') as f:
        pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)

    solution = simplex.solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=c0,
                                         maximization=maximization, varclass=varclass,
                                         pivotchoice=pivotchoice, cutpoint=cutpoint, callback=callback)

    sold = {'type':'dual', 'parameters':parameters, 'solution':solution}

    with open(solutionfilename, 'wb') as f:
        pickle.dump(sold, f, pickle.HIGHEST_PROTOCOL)

    return solution


def solveBirkhoffBothAndSave(primalproblemfilename=None, primalsolutionfilename=None,
                             dualproblemfilename=None, dualsolutionfilename=None,
                             lmax=2, k0 = 15, c = fractions.Fraction(17,10), *, varclass=fractions.Fraction,
                             pivotchoice=simplex.greedyStrategy, primalcutpoint=None, dualcutpoint=None,
                             callback=simplex.dummy, lset=None, mset=None,
                             kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                             fileprefix='', negativeBelly=False, negativeLeg=False, restrictionList=None,
                             jointLargeLeg=False, precomputedPartitionList=None,
                             Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer):
    """
    Creates, solves and saves both the primal and dual linear programs associated to the Birkhoff problem using
    the functions solveBirkhoffPrimalAndSave and solveBirkhoffDualAndSave. See these functions for meanings of
    parameters with the same name.
    - primalcutpoint is cutpoint used with solveBirkhoffPrimalAndSave.
    - dualcutpoint is cutpoint used with solveBirkhoffDualAndSave.
    - primalproblemfilename is problemfilename used with solveBirkhoffPrimalAndSave.
    - dualproblemfilename is problemfilename used with solveBirkhoffDualAndSave.
    - primalsolutionfilename is solutionfilename used with solveBirkhoffPrimalAndSave.
    - dualsolutionfilename is solutionfilename used with solveBirkhoffDualAndSave.

    Returns a tuple (xp, vp, xd, vd), where (xp, vp) is the return value of the first function and (xd, vd) is
    the return value of the second function.
    """
    
    (xp, vp) = solveBirkhoffPrimalAndSave(problemfilename=primalproblemfilename,
                                          solutionfilename=primalsolutionfilename,
                                          lmax=lmax, k0=k0, c=c, varclass=varclass,
                                          pivotchoice=pivotchoice, cutpoint=primalcutpoint, callback=callback,
                                          lset=lset, mset=mset, kostkaRounder=kostkaRounder,
                                          fileprefix=fileprefix,
                                          negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                          precomputedPartitionList=precomputedPartitionList)
    (xd, vd) = solveBirkhoffDualAndSave(problemfilename=dualproblemfilename,
                                        solutionfilename=dualsolutionfilename,
                                        lmax=lmax, k0=k0, c=c, varclass=varclass, 
                                        pivotchoice=pivotchoice, cutpoint=dualcutpoint, callback=callback,
                                        lset=lset, mset=mset, fileprefix=fileprefix,
                                        kostkaRounder=kostkaRounder, dualKostkaRounder=dualKostkaRounder,
                                        negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                        restrictionList=restrictionList,
                                        jointLargeLeg=jointLargeLeg,
                                        precomputedPartitionList=precomputedPartitionList,
                                        Yparametrizer=Yparametrizer,
                                        Wparametrizer=Wparametrizer)
    return (xp, vp, xd, vd)


def fragmentedPhaseFragBirkhoffDual(lmax, k0, c, *, varclass=fractions.Fraction,
                                    pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                                    lset=None, mset=None, kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                                    negativeBelly=False, negativeLeg=False,
                                    restrictionList=None, jointLargeLeg=False,
                                    save=True, fragproblemfilename=None, fragsolutionfilename=None,
                                    fileprefix='', precomputedPartitionList=None,
                                    Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer,
                                    parametersdict=None):
    """
    Creates and solves (and possibly saves) the fragmented phase of the fragment/defragment heuristic of the
    dual linear program associated to the Birkhoff problem using the functions createBirkhoffPrimal and
    simplex.solveNamedSimplex.

    If parameterdict is not None, uses parameterdict as entry 'parameters' of the saved problem and 
    fragmented solution, otherwise, creates an entry from parameters given to this function (parameterdict
    is mainly for use by other functions, e.g., solveFragBirkhoffDual, that wish to save different parameters
    into the 'parameters' entry).

    See solveFragBirkhoffDual for more details and meanings of parameters.

    Returns a pair (problem, sold), where problem is the problem in dictionary format and sold is the
    fragmented solution in dictionary format.
    """
    lsetstr = '' if lset==None else '_lset'
    msetstr = '' if mset==None else '_mset'
    if fragproblemfilename == None:
        fragproblemfilename = fileprefix + 'fragdualprob_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'
    if fragsolutionfilename == None:
        fragsolutionfilename = fileprefix + 'fragdualsol_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'

    if restrictionList == None:
        restrictionList = [(k, Partition.Partition([])) for k in range(1,k0+1)] + [(k, Partition.Partition([1,1])) for k in range(2,k0+1)] + [(1, Partition.Partition([l])) for l in range(1,lmax+1)]

    if parametersdict == None:
        parameters = {
            'lmax':lmax,
            'k0':k0,
            'c':c,
            'varclass':varclass,
            'pivotchoice':pivotchoice.__name__,
            'cutpoint':cutpoint,
            'callback':callback.__name__,
            'lset':lset,
            'mset':mset,
            'kostkaRounder':kostkaRounder.__name__,
            'dualKostkaRounder':dualKostkaRounder.__name__,
            'negativeBelly':negativeBelly,
            'negativeLeg':negativeLeg,
            'restrictionList':restrictionList,
            'jointLargeLeg':jointLargeLeg,
            'save':save,
            'fragproblemfilename':fragproblemfilename,
            'fragsolutionfilename':fragsolutionfilename,
            'fileprefix':fileprefix,
            'precomputedPartitionList':precomputedPartitionList,
            'Yparametrizer':Yparametrizer.__name__,
            'Wparametrizer':Wparametrizer.__name__,
            'parametersdict':parametersdict
        }
    else:
        parameters = parametersdict

    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffDual(lmax, k0, c, lset=lset, mset=mset,
                                                                     kostkaRounder=kostkaRounder,
                                                                     dualKostkaRounder=dualKostkaRounder,
                                                                     negativeBelly=negativeBelly,
                                                                     negativeLeg=negativeLeg,
                                                                     restrictionList=restrictionList,
                                                                     jointLargeLeg=jointLargeLeg,
                                                                     precomputedPartitionList=precomputedPartitionList,
                                                                     Yparametrizer=Yparametrizer,
                                                                     Wparametrizer=Wparametrizer)

    problem = {'c': c,
               'A': A,
               'b': b,
               'Aineq': Aineq,
               'vartype': vartype,
               'c0': c0,
               'maximization': maximization}

    d = {'type':'frag dual', 'parameters':parameters, 'problem':problem}

    if save:
        with open(fragproblemfilename, 'wb') as f:
            pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)

    solution = simplex.solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=c0,
                                         maximization=maximization, varclass=varclass,
                                         pivotchoice=pivotchoice, cutpoint=cutpoint, callback=callback)

    sold = {'type':'frag dual', 'parameters':parameters, 'solution':solution}

    if save:
        with open(fragsolutionfilename, 'wb') as f:
            pickle.dump(sold, f, pickle.HIGHEST_PROTOCOL)

    return (problem, sold)


def defragmentedPhaseFragBirkhoffDual(sold, *, save=True, defragsolutionfilename=None, fileprefix='',
                                      precomputedPartitionList=None,
                                      lightweight=1, verbose=False, verbosePeriod=100,
                                      num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50,
                                      addParameters=True, precision=21):
    """
    Runs (and possibly saves) the defragmented phase of the fragment/defragment heuristic of the
    dual linear program associated to the Birkhoff problem using the function inspectdualsol on the solution
    of the fragmented phase given in dictionary form by sold.

    The 'parameters' entry of the saved defragmented solution is the copied from the fragmented solution,
    except that entry 'restrictionList' is turned into 'original_restrictionList' and the new entry
    'restrictionList' is set to None. Further, if addParameters is True (default), then also adds parameters
    (other than sold) passed to this function to the 'parameters' entry (overwrites any entry that had the
    same name).

    See solveFragBirkhoffDual for more details and meanings of parameters.

    precision gives the precision of the printed objective value (if verbose is True).

    Returns defragmented solution in dictionary format.
    """
    parameters = sold['parameters']

    defragparameters = parameters.copy()        
    defragparameters['original_restrictionList'] = defragparameters['restrictionList']
    defragparameters['restrictionList'] = None

    if addParameters:
        defragParameters.update({
            'save':save,
            'defragsolutionfilename':defragsolutionfilename,
            'fileprefix':fileprefix,
            'precomputedPartitionList':precomputedPartitionList,
            'lightweight':lightweight,
            'verbose':verbose,
            'verbosePeriod':verbosePeriod,
            'num_workers':num_workers,
            'jobqueuemaxsize':jobqueuemaxsize,
            'resultqueuemaxsize':resultqueuemaxsize,
            'addParameters':addParameters
        })

    (x,v) = sold['solution']
    lmax = sold['parameters']['lmax']
    k0 = sold['parameters']['k0']
    c = sold['parameters']['c']

    mmax = 2*k0 + 2*lmax
    lset = list(range(0,lmax+1,2))
    mset = list(range(1,mmax+1))

    # Objective value computing
    obj = {}
    obj0 = 2
    for l in lset:
        obj['w[%03d]' % l] = -2 * tailbound(l, k0, c)

    for m in mset:
        obj['y[%03d]' % m] = - (c ** m - 1)
    
    (_, truevalue) = simplex.checkValue(c=obj, x=x, v=v, c0=obj0)

    defragsold = {'type':'defrag dual', 'parameters':defragparameters, 'solution':(x, truevalue)}

    if verbose:
        print('\nStarting defrag check.')

    insp = inspectdualsol(defragsold, lightweight=lightweight,
                          verbose=verbose, verbosePeriod=verbosePeriod,
                          precomputedPartitionList=precomputedPartitionList,
                          num_workers=num_workers,
                          jobqueuemaxsize=jobqueuemaxsize,
                          resultqueuemaxsize=resultqueuemaxsize)

    feasible = (len(insp[-1]) == 0)

    defragsold['feasible'] = feasible
    defragsold['insp'] = insp

    (_, truevalue) = simplex.checkValue(c=obj, x=x, v=v, c0=obj0)
    defragsold['solution'] = (x, truevalue)

    if save:
        with open(defragsolutionfilename, 'wb') as f:
            pickle.dump(defragsold, f, pickle.HIGHEST_PROTOCOL)

    return defragsold
    

def solveFragBirkhoffDual(lmax, k0, c, *, varclass=fractions.Fraction,
                          pivotchoice=simplex.greedyStrategy, cutpoint=None, callback=simplex.dummy,
                          lset=None, mset=None, kostkaRounder=dummy, dualKostkaRounder=dualDummy,
                          negativeBelly=False, negativeLeg=False, restrictionList=None, jointLargeLeg=False,
                          save=True, fragproblemfilename=None, fragsolutionfilename=None,
                          defragsolutionfilename=None, fileprefix='', precomputedPartitionList=None,
                          Yparametrizer=dummyParametrizer, Wparametrizer=dummyParametrizer,
                          dropPartitionListBeforeDefrag=False,
                          lightweight=1, verbose=False, verbosePeriod=100,
                          num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50,
                          runFragmentedPhaseOnly=False):
    """
    Creates and solves (and possibly saves) the fragment/defragment heuristic of the dual linear program
    associated to the Birkhoff problem using the functions createBirkhoffPrimal, simplex.solveNamedSimplex and
    inspectdualsol.
    
    The fragment/defragment heuristic is to solve a fragmented of the dual linear program that includes only a few
    restrictions and checking if the optimal solution of the fragmented is still a feasible (hence also optimal, 
    as long as no rounding was done in the fragmented version) solution of the full dual linear program (called
    defragmented dual linear program). Rounding is only performed in fragmented version of the problem.

    See simplex.solveNamedSimplex for meaning of parameters varclass, pivotchoice, cutpoint, callback.

    If verbose is True, prints messages for progress monitoring (this is also passed to inspectdualsol).
    See inspectdualsol for meaning of lightweight, verbosePeriod, num_workers, jobqueuemaxsize, resultqueuemaxsize.

    The other parameters are as follows.
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - kostkaRounder is a function to perform (dualization of the) primal rounding of Kostka coefficients.
    - dualKostkaRounder is a function to perform dual rounding of Kostka coefficients.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - restrictionList is a list (or generator) of which partitions should have their restrictions included. Each
    partition in this list is encoded as a tuple (k, belly), where k is its leg and belly is its belly. If this 
    parameter is None (default), then the function includes the following partitions in restrictionList.
    -- All hooks (b) of legs between 1 and k0 (inclusive).
    -- All b1_1 partitions of legs between 2 and k0 (inclusive).
    -- All partitions with leg 1 (hence belly of height 1) and belly size at most lmax.
    (Note that this default computation is different than the default computation of createBirkhoffDual and 
    partitionsRestrictionsBirkhoffDual.)
    - if jointLargeLeg is True, replaces all restrictions of corresponding to partitions (k, belly) with large
    leg (i.e., k > lmax) and same belly size s by a single restriction (k, s, lhs, ineq, rhs)
    that implies them all.
    - save: if True, then function also saves both fragmented problem, fragmented solution and
    defragmented solution on file. 
    - fragproblemfilename: name of file to save the fragmented problem (if None, a standardized name will be
    generated).
    - fragsolutionfilename: name of file to save the fragmented solution (if None, a standardized name will be
    generated).
    - defragsolutionfilename: name of file to save the defrag solution (if None, a standardized name will be
    generated).
    - fileprefix: string to be prepended to standardized names (if any is generated).
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.
    - Yparametrizer: if not dummyParametrizer, must be a fuction that parametrizes each y[m] variable according
    to the returned value of Yparametrizer(m). If Yparametrizer(m) == None, then y[m] is a free variable;
    else y[m] is made to satisfy
    y[m] = Yparametrizer(m) * Y + Yc,
    where Y and Yc are new variables (the same for all m).
    - Wparametrizer: if not dummyParametrizer, must be a fuction that parametrizes each w[l] variable according
    to the returned value of Wparametrizer(l). If Wparametrizer(l) == None, then w[l] is a free variable;
    else w[l] is made to satisfy
    w[l] = Wparametrizer(l) * W + Wc,
    where W and Wc are a new variables (the same for all l).
    - if runFragmentedPhaseOnly is True, then halts after finishing the fragmented phase instead; the returned
    defragmented solution dictionary is empty (see also resumeFragBirkhoffDual).

    WARNING: it is not a priori clear whether rounding (be it primal or dual) makes theoretical sense with the
    fragment/defragment heuristic.

    All file saves are made as dictionaries using the pickle module.

    Returns a triple (problem, sold, defragsold), where problem is the problem in dictionary format, sold is
    the fragmented solution in dictionary format and defragsold is the defragmented solution in dictionary
    format.
    """
    lsetstr = '' if lset==None else '_lset'
    msetstr = '' if mset==None else '_mset'
    if fragproblemfilename == None:
        fragproblemfilename = fileprefix + 'fragdualprob_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'
    if fragsolutionfilename == None:
        fragsolutionfilename = fileprefix + 'fragdualsol_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'
    if defragsolutionfilename == None:
        defragsolutionfilename = fileprefix + 'defragdualsol_' + str(lmax) + '_' + str(k0) + '_' + str(simplex.toDecimal(c)) + lsetstr + msetstr + '.pkl'

    if restrictionList == None:
        restrictionList = [(k, Partition.Partition([])) for k in range(1,k0+1)] + [(k, Partition.Partition([1,1])) for k in range(2,k0+1)] + [(1, Partition.Partition([l])) for l in range(1,lmax+1)]


    parameters = {
        'lmax':lmax,
        'k0':k0,
        'c':c,
        'varclass':varclass,
        'pivotchoice':pivotchoice.__name__,
        'cutpoint':cutpoint,
        'callback':callback.__name__,
        'lset':lset,
        'mset':mset,
        'kostkaRounder':kostkaRounder.__name__,
        'dualKostkaRounder':dualKostkaRounder.__name__,
        'negativeBelly':negativeBelly,
        'negativeLeg':negativeLeg,
        'restrictionList':restrictionList,
        'jointLargeLeg':jointLargeLeg,
        'save':save,
        'fragproblemfilename':fragproblemfilename,
        'fragsolutionfilename':fragsolutionfilename,
        'defragsolutionfilename':defragsolutionfilename,
        'fileprefix':fileprefix,
        'precomputedPartitionList':precomputedPartitionList,
        'Yparametrizer':Yparametrizer.__name__,
        'Wparametrizer':Wparametrizer.__name__,
        'lightweight':lightweight,
        'verbose':verbose,
        'verbosePeriod':verbosePeriod,
        'num_workers':num_workers,
        'jobqueuemaxsize':jobqueuemaxsize,
        'resultqueuemaxsize':resultqueuemaxsize,
        'runFragmentedPhaseOnly':runFragmentedPhaseOnly
    }

    (problem, sold) = fragmentedPhaseFragBirkhoffDual(lmax=lmax,
                                                      k0=k0,
                                                      c=c,
                                                      varclass=varclass,
                                                      pivotchoice=pivotchoice,
                                                      cutpoint=cutpoint,
                                                      callback=callback,
                                                      lset=lset,
                                                      mset=mset,
                                                      kostkaRounder=kostkaRounder,
                                                      dualKostkaRounder=dualKostkaRounder,
                                                      negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                                      restrictionList=restrictionList,
                                                      jointLargeLeg=jointLargeLeg,
                                                      save=save,
                                                      fragproblemfilename=fragproblemfilename,
                                                      fragsolutionfilename=fragsolutionfilename,
                                                      fileprefix=fileprefix,
                                                      precomputedPartitionList=precomputedPartitionList,
                                                      parametersdict=parameters,
                                                      Yparametrizer=Yparametrizer,
                                                      Wparametrizer=Wparametrizer)

    if runFragmentedPhaseOnly:
        return (problem, sold, dict())

    defragsold = defragmentedPhaseFragBirkhoffDual(sold,
                                                   save=save,
                                                   defragsolutionfilename=defragsolutionfilename,
                                                   fileprefix=fileprefix,
                                                   precomputedPartitionList=precomputedPartitionList,
                                                   lightweight=lightweight,
                                                   verbose=verbose,
                                                   verbosePeriod=verbosePeriod,
                                                   num_workers=num_workers,
                                                   jobqueuemaxsize=jobqueuemaxsize,
                                                   resultqueuemaxsize=resultqueuemaxsize,
                                                   addParameters=False)

    return (problem, sold, defragsold)


def resumeFragBirkhoffDual(sold):
    """
    Resumes computation of a call to solveFragBirkhoffDual with runFragmentedPhaseOnly=True from the halting
    point whose fragmented solution in dictionary format is sold.

    The following are different from running solveFragBirkhoffDual with runFragmentedPhaseOnly=False:
    - precomputedPartitionList is made to be None instead (this is would be equivalent to making
    dropPartitionListBeforeDefrag=True if Python was smart about freeing memory, but might actually be
    better than that if Python is not that smart).
    
    Returns defragmented solution in dictionary format.
    """
    save = sold['parameters']['save']
    defragsolutionfilename = sold['parameters']['defragsolutionfilename']
    fileprefix = sold['parameters']['fileprefix']
    precomputedPartitionList=None
    lightweight = sold['parameters']['lightweight']
    verbose = sold['parameters']['verbose']
    verbosePeriod = sold['parameters']['verbosePeriod']
    num_workers = sold['parameters']['num_workers']
    jobqueuemaxsize = sold['parameters']['jobqueuemaxsize']
    resultqueuemaxsize = sold['parameters']['resultqueuemaxsize']
    
    defragsold = defragmentedPhaseFragBirkhoffDual(sold,
                                                   save=save,
                                                   defragsolutionfilename=defragsolutionfilename,
                                                   fileprefix=fileprefix,
                                                   precomputedPartitionList=precomputedPartitionList,
                                                   lightweight=lightweight,
                                                   verbose=verbose,
                                                   verbosePeriod=verbosePeriod,
                                                   num_workers=num_workers,
                                                   jobqueuemaxsize=jobqueuemaxsize,
                                                   resultqueuemaxsize=resultqueuemaxsize,
                                                   addParameters=False)

    return defragsold


#######################
#  File manipulators  #
#######################

def load(filename):
    """
    Loads and returns pickle object from file filename.
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)

def loadall(foldername, extension='.pkl', *, sortkey=lambda s : (s['type'], s['parameters']['lmax'])):
    """
    Loads all pickle objects from files with extension extension from folder foldername (without trailing slash)
    into a list ret, sorts ret using sortkey and returns ret.
    """
    ret = []
    for f in glob.glob(foldername + '/*' + extension):
        ret.append(load(f))
    ret.sort(key=sortkey)
    return ret

#######################
#  Solution printers  #
#######################

def printDualNormalizedVariables(x, c):
    """
    Prints an attempt at a normalization of variables x from a dual solution of a problem with
    pseudorandomness parameter c.
    """
    l = [(s, v) for (s,v) in x.items() if v != 0]
    l.sort()
    for (s,v) in l:
        mult = 1
        if s[0] == 'y':
            mult = c ** int(s[2:-1])
        print(s + ' = %s' % str(simplex.toDecimal(v * mult)))

def printsolution(x, v, precision=21, padding=10):
    """
    Prints solution (x,v) in Decimal format with precision precision and padding padding.
    """
    l = [(name, value) for (name, value) in x.items() if value != 0]
    l.sort()
    for (name, value) in l:
        print('%*s = %s' % (padding, name, str(simplex.toDecimal(value, precision))))
    print('%*s = %s' % (padding, 'Value', str(simplex.toDecimal(v, precision))))

def printrawsolution(x, v, padding=10):
    """
    Prints solution (x,v) in raw format with padding padding.
    """
    l = [(name, value) for (name, value) in x.items() if value != 0]
    l.sort()
    for (name, value) in l:
        print('%*s = ' % (padding, name) + str(value))
    print('%*s = ' % (padding, 'Value') + str(v))

def printsol(sol, precision=21, padding=10):
    """
    Prints dictionary formatted (see solveBirkhoffPrimalAndSave) solution sol in Decimal format with
    precision precision and padding padding.
    """
    soltype = sol.get('type', None)
    lmax = k0 = c = None
    parameters = sol.get('parameters', None)
    feasible = sol.get('feasible', None)
    if parameters != None:
        lmax = parameters.get('lmax', None)
        k0 = parameters.get('k0', None)
        c = parameters.get('c', None)

    if soltype != None:
        print('%*s = %s' % (padding, 'type', soltype))
    if feasible != None:
        print('%*s = %s' % (padding, 'feasible', feasible))
    if lmax != None:
        print('%*s = %d' % (padding, 'lmax', lmax))
    if k0 != None:
        print('%*s = %d' % (padding, 'k0', k0))
    if c != None:
        print('%*s = %s' % (padding, 'c', str(simplex.toDecimal(c, precision))))

    return printsolution(*sol['solution'], precision=precision, padding=padding)

def printrawsol(sol, padding=10):
    """
    Prints dictionary formatted (see solveBirkhoffPrimalAndSave) solution sol in raw format with padding padding.
    """
    soltype = sol.get('type', None)
    lmax = k0 = c = None
    parameters = sol.get('parameters', None)
    feasible = sol.get('feasible', None)
    if parameters != None:
        lmax = parameters.get('lmax', None)
        k0 = parameters.get('k0', None)
        c = parameters.get('c', None)

    if soltype != None:
        print('%*s = %s' % (padding, 'type', soltype))
    if feasible != None:
        print('%*s = %s' % (padding, 'feasible', feasible))
    if lmax != None:
        print('%*s = %d' % (padding, 'lmax', lmax))
    if k0 != None:
        print('%*s = %d' % (padding, 'k0', k0))
    if c != None:
        print('%*s = %s' % (padding, 'c', str(c)))

    return printrawsolution(*sol['solution'], padding=padding)

def printallsol(sols, precision=21, padding=10, separator='\n'):
    """
    Prints all solutions in sols using printsol with parameters precision and padding.
    Separates solutions using separator.
    """
    for s in sols:
        printsol(s, precision=precision, padding=padding)
        print(separator, end='')

def printallrawsol(sols, padding=10, separator='\n'):
    """
    Prints all solutions in sols using printrawsol with parameters precision and padding.
    Separates solutions using separator.
    """
    for s in sols:
        printrawsol(s, padding=padding)
        print(separator, end='')

#####################
#  Solution fixers  #
#####################

def dontFix(sol, k, belly, rest, ineq, rhs):
    """
    A solution fixer that does not do anything.
    """
    return False

def singleyFix(sol, k, belly, rest, ineq, rhs):
    """
    A solution fixer that increases the best y the least amount possible.
    """
    if ineq != '>=' or rhs != 0:
        return False

    c = sol['parameters']['c']
    (x,v) = sol['parameters']['solution']
    bestvar = None
    bestpayoff = 0
    val = 0
    for (var, coeff) in rest.items():
        if var[0] == 'y':
            m = int(var[1:-1])
            payoff = coeff / (c**m - 1)
            if payoff > bestpayoff:
                bestpayoff = payoff
                bestvar = var
        xvar = x.get(var, 0)
        if xvar != 0 and coeff != 0:
            val += coeff * xvar
    if val < 0:
        x[bestvar] -= val / rest[bestvar]
        
    return False


#########################
#  Solution inspectors  #
#########################

def checkprobsol(prob, sol):
    """
    Checks feasibility and correctness of objective value of solution sol for problem prob, both in dictionary
    format (see solveBirkhoffPrimalAndSave) using functions simplex.checkFeasibility and simplex.checkValue.

    Returns tuple (correct, slack, truevalue), where correct is True if and only if sol is both feasible and
    has correct value, slack is as returned by simplex.checkFeasibility and truevalue is as returned by
    simplex.checkValue.
    """

    problem = prob['problem']
    c = problem['c']
    A = problem['A']
    b = problem['b']
    c0 = problem['c0']
    Aineq = problem['Aineq']
    vartype = problem['vartype']

    (x,v) = sol['solution']

    (feasible, slack) = simplex.checkFeasibility(A,b, Aineq, vartype, x)
    (valcorrect, truevalue) = simplex.checkValue(c, x, v, c0=c0)

    return (feasible and valcorrect, slack, truevalue)

def checksol(sol):
    """
    Checks feasibility and correctness of objective value of solution sol in dictionary
    format (see solveBirkhoffPrimalAndSave) by first generating the problem again using createBirkhoffPrimal,
    then using functions simplex.checkFeasibility and simplex.checkValue.

    WARNING: if rounding was performed for the original problem to which sol is associated, the problem 
    generated by createBirkhoffPrimal will not match the original one.

    Returns tuple (correct, slack, truevalue), where correct is True if and only if sol is both feasible and
    has correct value, slack is as returned by simplex.checkFeasibility and truevalue is as returned by
    simplex.checkValue.
    """

    parameters = sol['parameters']
    lmax = parameters['lmax']
    k0 = parameters['k0']
    c = parameters['c']
    lset = parameters['lset']
    mset = parameters['mset']
    negativeBelly = parameters['negativeBelly']
    negativeLeg = parameters['negativeLeg']

    (x,v) = sol['solution']

    (c, A, b, Aineq, vartype, c0, maximization) = createBirkhoffPrimal(lmax, k0, c,
                                                                       lset=lset, mset=mset,
                                                                       kostkaRounder=dummy,
                                                                       negativeBelly=negativeBelly,
                                                                       negativeLeg=negativeLeg)

    (feasible, slack) = simplex.checkFeasibility(A,b, Aineq, vartype, x)
    (valcorrect, truevalue) = simplex.checkValue(c, x, v, c0=c0)

    return (feasible and valcorrect, slack, truevalue)

def checkallsol(sols, *, verbose=True):
    """
    Checks feasibility and correctness of all solutions in sols using checksol.

    If verbose is True, prints messages for progress monitoring.

    Returns list with respective return values.
    """
    ret = []
    for i in range(len(sols)):
        sol = sols[i]
        if verbose:
            print('Checking %2d of %2d:' % (i+1, len(sols)))
            soltype = sol.get('type', None)
            lmax = k0 = c = None
            parameters = sol.get('parameters', None)
            if parameters != None:
                lmax = parameters.get('lmax', None)
                k0 = parameters.get('k0', None)
                c = parameters.get('c', None)

            if soltype != None:
                print('%s = %s' % ('type', soltype))
            if lmax != None:
                print('%s = %d' % ('lmax', lmax))
            if k0 != None:
                print('%s = %d' % ('k0', k0))
            if c != None:
                print('%s = %s' % ('c', str(c)))

        ret.append(checksol(sol))

        if verbose:
            print(('Correct!' if ret[-1][0] else 'Not correct!'), end='\n\n')

    return ret

class FakeLock:
    def __enter__(self):
        self.acquire()
        return self
    def __exit__(self, tp, val, traceback):
        self.release()
    def acquire(self, block=True, timeout=None):
        return True
    def release(self):
        pass


def inspectdualsol(sol, *, lightweight=0, verbose=False, verbosePeriod=100, precision=21, restrictionList=None,
                   precomputedPartitionList=None, solutionFixer=dontFix,
                   num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50):
    """
    Inspects which dual partition restrictions are satisfied by weight alone, which have slack and which are tight.
    Returns a tuple (veryeasy, easy, medium, hard, notsat); if lightweight is 0, each of these is a
    dictionary whose entries are (name, k):(wval, tval), where (name, k) represents a partition with belly b such
    that b.bname() is name and leg k, wval is the value of the weight part of the restriction and tval is the
    total value of the restriction; if lightweight is 1, each of these is a CompactSet
    representing a set containing only the (name, k) part of the lightweight=False case; if lightweight is 2, the
    entries veryeasy and medium are made to be the empty set (since these are typically very large).

    Furthermore, we have
    - veryeasy contains the partitions that are satisfied by weight alone and have slackness.
    - easy contains the partitions whose restrictions are satisfied by weight alone and are tight.
    - medium contains the partitions whose restrictions are satisfied by using weight and Young but have slackness.
    - hard contains the partitions whose restrictions are satisfied by using weight and Young and are tight.
    - notsat contains the partitions whose restrictions are not satisfied at all.

    The other parameters of this function are as follows.
    - verbose: if True, prints messages for progress monitoring.
    - verbosePeriod gives the period of monitoring messages.
    - if restrictionList is not None, checks only the restrictions in restrictionList; otherwise, reads
    restrictionList from sol['parameters']['restrictionList'] and uses that instead; if this entry does not exist
    or is None, then the function includes all partitions of belly with size at most lmax and leg
    at most k0.
    - precision is the precision of the printed objective value (if verbose is True).
    - precomputedPartitionList is assumed to be one of the following.
    - solutionFixer is a function to be called whenever a restriction that was not satisfied is found to attempt
    to fix the solution (this means that sol can be changed). The parameters passed to solutionFixer are
    -- sol: dictionary of the solution
    -- k: leg of restriction that is not satisfied
    -- belly: belly of restriction that is not satisfied
    -- rest: left-hand side of restriction that is not satisfied
    -- ineq: inequality type of restriction that is not satisfied
    -- rhs: right-hand side of restriction that is not satisfied
    The return value of solutionFixer should be False if the solution fix is conservative in the sense that it
    cannot have made restrictions that were previously satisfied not satisfied; otherwise it should be True (in
    this case inspectdualsol will restart the solution check). If verbose is True, after the fix is performed the
    objective value is recomputed and printed again. If such fix is performed, the returned value may not be
    accurate. The function solutionFixer should change entries in sol rather than reassigning them new containers.
    The default solutionFixer function dontFix does not do anything.
    -- A list of length at least lmax+1 such that precomputedPartitionList[l] is the list of all
    partitions of size l.
    -- None, in which case such list will be computed by this function.
    - num_workers: if 1, executes normally, if a number larger than 1, executes using the multiprocessing module
    with num_workers workers (since one worker is allocated to generate the restrictions, this should be at least
    3 to make a performance difference).
    - jobqueuemaxsize: maximum size of queue for producer (in multiprocessing mode)
    - resultqueuemaxsize: maximum size of queue for collector (in multiprocessing mode)

    WARNING: if rounding was performed for the original problem to which sol is associated, the restrictions
    generated by this function will not match the original ones.
    """

    parameters = sol['parameters']
    (x, v) = sol['solution']

    lmax = parameters['lmax']
    k0 = parameters['k0']
    c = parameters['c']
    lset = parameters.get('lset', None)
    mset = parameters.get('mset', None)
    negativeBelly = parameters.get('negativeBelly', False)
    negativeLeg = parameters.get('negativeLeg', False)
    if restrictionList == None:
        restrictionList = parameters.get('restrictionList', None)
    jointLargeLeg = parameters.get('jointLargeLeg', False)

    mmax = 2*k0 + 2*lmax

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    if verbose:
        # Objective value computing
        obj = {}
        obj0 = 2
        for l in lset:
            obj['w[%03d]' % l] = -2 * tailbound(l, k0, c)

        for m in mset:
            obj['y[%03d]' % m] = - (c ** m - 1)

        (_, truevalue) = simplex.checkValue(c=obj, x=x, v=v, c0=obj0)
        print('Objective value: %s ' % str(simplex.toDecimal(truevalue, precision=precision)))


    if num_workers == 1: # Run without multiprocessing
        partitions = None
        if verbose:
            print('Computing approximate total number of checks required. . .')
            if restrictionList == None:
                partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList

                totalChecks = 0
                for i in range(len(partitions)):
                    print('\rSteps concluded: %3d of %d' % (i, len(partitions)), end='')
                    if jointLargeLeg:
                        totalChecks += sum(((lmax + 1 - height) * Partition.countPartitions(i, width=height) for height in range(i+1))) + (k0 - lmax) * (lmax + 1)
                    else:
                        totalChecks += sum(((k0 + 1 - height) * Partition.countPartitions(i, width=height) for height in range(i+1)))

                print('\rSteps concluded: %3d of %d' % (len(partitions), len(partitions)))                    
            else:
                totalChecks = len(restrictionList)

        while True:
            repeat = False

            partRests = partitionRestrictionsBirkhoffDual(lmax=lmax, k0=k0, lset=lset, mset=mset,
                                                          negativeBelly=negativeBelly, negativeLeg=negativeLeg,
                                                          restrictionList=restrictionList,
                                                          jointLargeLeg=jointLargeLeg,
                                                          precomputedPartitionList=partitions)

            ret = (CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet()) if lightweight else (dict(), dict(), dict(), dict(), dict())

            if verbose:
                cnt = [0] * 5
                vcnt = 0
                print('\n', end='')
            for (k, belly, rest, ineq, rhs) in partRests:
                tval = wval = 0
                for (key, restval) in rest.items():
                    xkey = x.get(key, 0)
                    if xkey != 0 and restval != 0:
                        if key[0] == 'w':
                            wval += restval * xkey
                        else:
                            tval += restval * xkey
                tval += wval

                ind = 2 * (tval < 0) + 2 * (wval < 0) + (tval == 0)
                if lightweight:
                    if lightweight == 1 or (ind != 0 and ind != 2):
                        ret[ind].add((belly.bname(), k))
                else:
                    ret[ind][(belly.bname(), k)] = (wval, tval)
                if ind == 4:
                    if verbose:
                        print('\n%s is not satisfied' % str((belly.bname(), k)))
                    recheck = solutionFixer(sol, k, belly, rest, ineq, rhs)
                    if verbose:
                        (_, truevalue) = simplex.checkValue(c=obj, x=x, v=v, c0=obj0)
                        print('Fixed objective value: %s ' % str(simplex.toDecimal(truevalue, precision=precision)))
                    if reckeck:
                        repeat = True
                        if verbose:
                            print('Recheck required.\n')
                        break
                if verbose:
                    cnt[ind] += 1
                    vcnt += 1
                    if vcnt % verbosePeriod == 0:
                        print('\rAp.total: %12d v.easy: %12d, easy: %3d, medium: %12d, hard: %3d, notsat: %3d' % ((totalChecks,) + tuple(cnt)), end='')
            if repeat == False:
                break

        if verbose:
            print('\rAp.total: %12d v.easy: %12d, easy: %3d, medium: %12d, hard: %3d, notsat: %3d' % ((totalChecks,) + tuple(cnt)))
            print(('Feasible!' if cnt[-1] == 0 else 'Not feasible!'))

        if lightweight:
            for cs in ret:
                cs.reorganize()
            ret = (CompactSet.CompactSet(ret[0]),
                   CompactSet.CompactSet(ret[1]),
                   CompactSet.CompactSet(ret[2]),
                   CompactSet.CompactSet(ret[3]),
                   CompactSet.CompactSet(ret[4]))

        return ret

    # run using multiprocessing

    if verbose:
        partitions = [Partition.enumeratePartitions(l) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList
        # We use generators here to save memory (so they need to be reinitialized)
        print('Computing approximate total number of checks required. . .')
        totalChecks = 0
        if restrictionList == None:
            for i in range(len(partitions)):
                print('\rSteps concluded: %3d of %d' % (i, len(partitions)), end='')
                if jointLargeLeg:
                    totalChecks += sum(((lmax + 1 - height) * Partition.countPartitions(i, width=height) for height in range(i+1))) + (k0 - lmax) * (lmax + 1)
                else:
                    totalChecks += sum(((k0 + 1 - height) * Partition.countPartitions(i, width=height) for height in range(i+1)))
            print('\rSteps concluded: %3d of %d' % (len(partitions), len(partitions)))
        else:
            totalChecks = len(restrictionList)

    #  Multiprocessing functions  #

    def produce(num_consumers, jobqueue, restrictionList, precomputedPartitionList):
        partitions = [Partition.enumeratePartitions(l) for l in range(lmax+1)] if precomputedPartitionList == None else precomputedPartitionList
        # We use generators here to save memory (so they need to be reinitialized)

        if restrictionList == None:
            restrictionGen = ((k, belly) for bellysize in range(0, lmax+1) for belly in partitions[bellysize] for k in range(max(belly.height(), 1), k0+1))
        else:
            restrictionGen = (rest for rest in restrictionList)

        largeLegPairs = set()

        for (k, belly) in restrictionGen:
            if k > lmax and jointLargeLeg:
                largeLegPairs.add((k, Partition.JointPartition(belly.size())))
                continue

            jobqueue.put((k, belly))


        for (k, jp) in largeLegPairs: # if jointLargeLeg == False, this is empty
            if jp.size() in lset and (k % 2 == 1 or negativeLeg):
                # Adds only if Parseval would give some negative coefficient.
                # (If all coefficients are non-negative, then the restriction is trivially satisfied.)

                jobqueue.put((k, jp))

        for i in range(num_consumers):
            jobqueue.put(None) # A signal for consumers that the task is complete

    def consume(jobqueue, resultqueue, lock, readxflag, sharedx):
        while True:
            item = jobqueue.get()
            with lock:
                if readxflag.is_set(): 
                    x.clear()
                    x.update(sharedx)
                    readxflag.clear()

                if item == None:
                    break

                (k, belly) = item
                if type(belly) == Partition.Partition:
                    (_, _, rest, ineq, rhs) = generateSinglePartitionRestriction(k=k, belly=belly,
                                                                                 lset=lset, mset=mset,
                                                                                 negativeBelly=negativeBelly,
                                                                                 negativeLeg=negativeLeg)
                else: # type(belly) == Partition.JointPartition
                    (_, _, rest, ineq, rhs) = generateSingleLargeLegPartitionRestriction(k=k, belly=belly,
                                                                                         lset=lset, mset=mset,
                                                                                         negativeLeg=negativeLeg)

                tval = wval = 0
                for (key, restval) in rest.items():
                    xkey = x.get(key, 0)
                    if xkey != 0 and restval != 0:
                        if key[0] == 'w':
                            wval += restval * xkey
                        else:
                            tval += restval * xkey
                tval += wval

                ind = 2 * (tval < 0) + 2 * (wval < 0) + (tval == 0)

            # We release the lock before putting the result in the queue to avoid a deadlock
            if lightweight:
                resultqueue.put((ind, (belly, k)))
            else:
                resultqueue.put((ind, (belly, k), (wval, tval)))

        resultqueue.put(None) # A signal for the collector that the task is complete


    def collect(num_consumers, jobqueue, resultqueue, conlocks, conreadxflags, sharedx):
        signals = 0
        ret = (CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet(), CompactSet.SemiCompactSet()) if lightweight else (dict(), dict(), dict(), dict(), dict())

        if verbose:
            cnt = [0] * 5
            vcnt = 0
        
        while signals < num_consumers:
            item = resultqueue.get()
            if item == None:
                signals += 1
                continue
            if lightweight:
                (ind, (belly, k)) = item
                if lightweight == 1 or (ind != 0 and ind != 2):
                    ret[ind].add((belly.bname(), k))
            else:
                (ind, (belly, k), val) = item
                ret[ind][(belly.bname(), k)] = val
            if verbose:
                cnt[ind] += 1
                vcnt += 1
                if vcnt % verbosePeriod == 0:
                    print('\rAp.total: %12d v.easy: %12d, easy: %3d, medium: %12d, hard: %3d, notsat: %3d' % ((totalChecks,) + tuple(cnt)), end='')

            if ind == 4:
                if verbose:
                    print('%s is not satisfied' % str((belly.bname(), k)))

                if solutionFixer != dontFix:
                    for lock in conlocks:
                        lock.acquire() # Stops all consumers
                    for flag in conreadxflags:
                        flag.set() # Signals the consumers that x must be read from sharedx

                    # Regenerating restriction for fixer
                    if type(belly) == Partition.Partition:
                        (_, _, rest, ineq, rhs) = generateSinglePartitionRestriction(k=k, belly=belly,
                                                                                     lset=lset, mset=mset,
                                                                                     negativeBelly=negativeBelly,
                                                                                     negativeLeg=negativeLeg)
                    else: # type(belly) == Partition.JointPartition
                        (_, _, rest, ineq, rhs) = generateSingleLargeLegPartitionRestriction(k=k, belly=belly,
                                                                                             lset=lset, mset=mset,
                                                                                             negativeLeg=negativeLeg)

                    recheck = solutionFixer(sol, k, belly, rest, ineq, rhs)
                    if verbose:
                        (_, truevalue) = simplex.checkValue(c=obj, x=x, v=v, c0=obj0)
                        print('Fixed objective value: %s ' % str(simplex.toDecimal(truevalue, precision=precision)))

                    sharedx.clear()
                    sharedx.update(x)

                    if recheck:
                        if verbose:
                            print('Recheck required.\n')
                        return None

                    for lock in conlocks:
                        lock.release() # Resumes all consumers

        if verbose:
            print('\rAp.total: %12d v.easy: %12d, easy: %3d, medium: %12d, hard: %3d, notsat: %3d' % ((totalChecks,) + tuple(cnt)))
            print(('Feasible!' if cnt[-1] == 0 else 'Not feasible!'))

        if lightweight:
            for cs in ret:
                cs.reorganize()
            ret = (CompactSet.CompactSet(ret[0]),
                   CompactSet.CompactSet(ret[1]),
                   CompactSet.CompactSet(ret[2]),
                   CompactSet.CompactSet(ret[3]),
                   CompactSet.CompactSet(ret[4]))

        return ret

    #  End of multiprocessing functions  #

    num_consumers = num_workers - 1

    while True:
        try:
            jobqueue = multiprocessing.Queue(maxsize=jobqueuemaxsize)
            resultqueue = multiprocessing.Queue(maxsize=resultqueuemaxsize)

            producer = multiprocessing.Process(target=produce, args=(num_consumers,
                                                                     jobqueue,
                                                                     restrictionList,
                                                                     precomputedPartitionList))
            producer.start()

            conlocks = [(multiprocessing.Lock() if solutionFixer!=dontFix else FakeLock())
                        for i in range(num_consumers)]
            conreadxflags = [multiprocessing.Event() for i in range(num_consumers)]
            manager = multiprocessing.Manager()
            sharedx = manager.dict()

            consumers = [multiprocessing.Process(target=consume, args=(jobqueue,
                                                                       resultqueue,
                                                                       conlocks[i],
                                                                       conreadxflags[i],
                                                                       sharedx))
                         for i in range(num_consumers)]

            for c in consumers:
                c.start()

            ret = collect(num_consumers,
                          jobqueue,
                          resultqueue,
                          conlocks,
                          conreadxflags,
                          sharedx)
            if ret == None:
                for c in consumers:
                    c.terminate()
                producer.terminate()
            else:
                break
        finally:
            for c in consumers:
                c.join()
            producer.join()

    return ret


def inspectalldualsol(sols, *, verbose=True, verbosePeriod=100, precomputedPartitionList=None,
                      num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50):
    """
    Inspects all dual solutions in sols using inspectdualsol.

    If verbose is True, prints messages for progress monitoring (this is also passed to inspectdualsol).
    verbosePeriod, precomputedPartitionList, num_workers, jobqueuemaxsize, resultqueuemaxsize are passed
    to inspectdualsol.

    Returns list with respective return values.
    """
    ret = []
    dsols = [sol for sol in sols if sol.get('type', '')[-4:] == 'dual']
    for i in range(len(dsols)):
        sol = dsols[i]
        if verbose:
            print('\nInspecting %2d of %2d:' % (i+1, len(dsols)))
            soltype = sol.get('type', None)
            lmax = k0 = c = None
            parameters = sol.get('parameters', None)
            if parameters != None:
                lmax = parameters.get('lmax', None)
                k0 = parameters.get('k0', None)
                c = parameters.get('c', None)

            if soltype != None:
                print('%s = %s' % ('type', soltype))
            if lmax != None:
                print('%s = %d' % ('lmax', lmax))
            if k0 != None:
                print('%s = %d' % ('k0', k0))
            if c != None:
                print('%s = %s' % ('c', str(c)))

        ret.append(inspectdualsol(sol, verbose=verbose,
                                  verbosePeriod=verbosePeriod,
                                  precomputedPartitionList=precomputedPartitionList,
                                  num_workers=num_workers,
                                  jobqueuemaxsize=jobqueuemaxsize,
                                  resultqueuemaxsize=resultqueuemaxsize))

    return ret

def computeAverageFractionSize(sol):
    """
    Returns the average approximate size of the non-zero fractions appearing in solution sol using
    simplex.sizeOfFraction.
    """
    (x,v) = sol['solution']

    num = simplex.sizeOfFraction(v)
    den = 1
    for (key, value) in x.items():
        if value != 0:
            num += simplex.sizeOfFraction(value)
            den += 1

    return num / den

def computeMaxFractionSize(sol):
    """
    Returns the maximum approximate size of the fractions appearing in solution sol using
    simplex.sizeOfFraction.
    """
    (x,v) = sol['solution']

    maxsize = simplex.sizeOfFraction(v)
    for (key, value) in x.items():
        size = simplex.sizeOfFraction(value)
        if size > maxsize:
            maxsize = size

    return maxsize

def computeTotalFractionSize(sol):
    """
    Returns the maximum approximate size of the fractions appearing in solution sol using
    simplex.sizeOfFraction.
    """
    (x,v) = sol['solution']

    total = simplex.sizeOfFraction(v)
    for (key, value) in x.items():
        total += simplex.sizeOfFraction(value)

    return total


#############################
#  Solution simplification  #
#############################

def generateDecimalSimplifier(*, precision=30, rounding=decimal.ROUND_CEILING, Emin=-999999999, Emax=999999999):
    """
    Generates a decimal simplifier with precision precision, rounding rounding, Emin Emin and Emax Emax
    (see decimal.Context).

    The simplifier argument x is assumed to represent the number x.numerator / x.denominator.
    """
    def simplifier(x):
        with decimal.localcontext(decimal.Context(prec=precision, rounding=rounding, Emin=Emin, Emax=Emax)):
            decx = decimal.Decimal(x.numerator) / decimal.Decimal(x.denominator)
        return type(x)(decx)
    return simplifier

def simplifySolution(sol, simplifier):
    """
    Simplifies solution sol in dictionary format using function simplifier.
    Returns simplified solution in dictionary format.

    No guarantees are given as to feasibility or objective value, except that weights on dual solutions are
    forced to remain adding to 1.

    Auxiliary parameters Y, Yc, W or Wc for Yparametrizers and Wparametrizers are removed from
    the solution.
    """
    simpsol = sol.copy()
    (x,v) = sol['solution']
    simpv = simplifier(v)
    simpx = {key:simplifier(val) for (key, val) in x.items()
             if key != 'Y' and key != 'Yc' and key != 'W' and key != 'Wc'}
    simpsol['solution'] = (simpx, simpv)

    if sol['type'] == 'primal':
        return simpsol

    if sol['type'] == 'dual' or sol['type'] == 'frag dual' or sol['type'] == 'defrag dual':
        maxval = 0
        maxkey = None
        weightsum = 0
        for (key, val) in simpx.items():
            if key[0] == 'w':
                weightsum += val
                if val > maxval:
                    maxval = val
                    maxkey = key
        simpx[maxkey] += 1 - weightsum
        return simpsol

    raise Exception('Unknown solution type: %s' % sol['type'])


def simplifyFile(filename, simplifier, savefileprefix='', *, savefilename=None, copyfiletimes=True):
    """
    Reads a dictionary formatted solution (see solveFragBirkhoffDual) from file filename, applies
    simplifySolution to it and saves it in file savefilename.

    If savefilename is None (default), then uses savefileprefix + filename as savefilename.

    If copyfiletimes is True, the access time and modified time of file filename will be copied to file 
    savefilename.
    """
    if savefilename == None:
        savefilename = savefileprefix + filename

    sol = load(filename)
    if copyfiletimes:
        filestat = os.stat(filename)
    
    simpsol = simplifySolution(sol, simplifier)
    with open(savefilename, 'wb') as f:
        pickle.dump(simpsol, f, pickle.HIGHEST_PROTOCOL)

    if copyfiletimes:
        os.utime(savefilename, (filestat.st_atime, filestat.st_mtime))


def simplifyAllFiles(foldername, simplifier, extension='.pkl', *, copyfiletimes=True, verbose=True):
    """
    Applies simplifyFile to all files from folder foldername (without trailing slash) with extension extension.
    The parameter savefilename from simplifyFile used is None and the parameter savefileprefix is ''.
    The parameters simplifier and copyfiletimes are passed on to simplifyFile.

    If verbose is True, prints messages for monitoring.
    """
    for f in glob.glob(foldername + '/*' + extension):
        if verbose:
            print("Processing file %s. . ." % f)
        simplifyFile(filename=f, simplifier=simplifier, savefileprefix='', copyfiletimes=copyfiletimes)

def simplifyFixSimplifySolution(sol, simplifier, *, verbose=True, verbosePeriod=1, precision=21,
                                restrictionList=None, precomputedPartitionList=None, solutionFixer=singleyFix,
                                num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50):
    """
    Applies simplifySolution, then inspectdualsol, then simplifySolution to a solution.
    The purpose of inspectdualsol is to fix the solution.

    See simplifySolution and inspectdualsol for meanings of parameters. The only change is on the parameter
    restrictionList, which when None is made to be the following list of partitions.
    -- All hooks (b) of legs between 1 and k0 (inclusive).
    -- All b1_1 partitions of legs between 2 and k0 (inclusive).
    -- All partitions with leg 1 (hence belly of height 1) and belly size at most lmax.
    
    Returns final solution.
    """
    lmax = sol['parameters']['lmax']
    k0 = sol['parameters']['k0']
    
    if restrictionList == None:
        restrictionList = [(k, Partition.Partition([])) for k in range(1,k0+1)] + [(k, Partition.Partition([1,1])) for k in range(2,k0+1)] + [(1, Partition.Partition([l])) for l in range(1,lmax+1)]

    
    newsol = simplifySolution(sol, simplifier)

    inspectdualsol(newsol, lightweight=2, verbose=verbose, verbosePeriod=verbosePeriod, precision=precision,
                   restrictionList=restrictionList, precomputedPartitionList=precomputedPartitionList,
                   solutionFixer=solutionFixer, num_workers=num_workers, jobqueuemaxsize=jobqueuemaxsize,
                   resultqueuemaxsize=resultqueuemaxsize)

    return simplifySolution(newsol, simplifier)

def simplifyFixSimplifyFile(filename, simplifier, savefileprefix='', *, savefilename=None, copyfiletimes=True,
                            verbose=True, verbosePeriod=1, precision=21,
                            restrictionList=None, precomputedPartitionList=None, solutionFixer=singleyFix,
                            num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50):
    """
    Reads a dictionary formatted solution (see solveFragBirkhoffDual) from file filename, applies
    simplifyFixSimplifySolution to it and saves it in file savefilename.

    If savefilename is None (default), then uses savefileprefix + filename as savefilename.

    If copyfiletimes is True, the access time and modified time of file filename will be copied to file 
    savefilename.

    See simplifyFixSimplifySolution for meanings of other parameters.
    """
    if savefilename == None:
        savefilename = savefileprefix + filename

    sol = load(filename)
    if copyfiletimes:
        filestat = os.stat(filename)
    
    simpsol = simplifyFixSimplifySolution(sol, simplifier,
                                          verbose=verbose,
                                          verbosePeriod=verbosePeriod,
                                          precision=precision,
                                          restrictionList=restrictionList,
                                          precomputedPartitionList=precomputedPartitionList,
                                          solutionFixer=solutionFixer,
                                          num_workers=num_workers,
                                          jobqueuemaxsize=jobqueuemaxsize,
                                          resultqueuemaxsize=resultqueuemaxsize)
    with open(savefilename, 'wb') as f:
        pickle.dump(simpsol, f, pickle.HIGHEST_PROTOCOL)

    if copyfiletimes:
        os.utime(savefilename, (filestat.st_atime, filestat.st_mtime))


def simplifyFixSimplifyAllFiles(foldername, simplifier, extension='.pkl', *, copyfiletimes=True, verbose=True,
                                verbosePeriod=1, precision=21,
                                restrictionList=None, precomputedPartitionList=None, solutionFixer=singleyFix,
                                num_workers=1, jobqueuemaxsize=50, resultqueuemaxsize=50):
    """
    Applies simplifyFixSimplifyFile to all files from folder foldername (without trailing slash) with
    extension extension.
    The parameter savefilename from simplifyFile used is None and the parameter savefileprefix is ''.
    The other parameters are passed on to simplifyFile.

    If verbose is True, prints messages for monitoring.
    """
    for f in glob.glob(foldername + '/*' + extension):
        if verbose:
            print("Processing file %s. . ." % f)
        simplifyFixSimplifyFile(filename=f,
                                simplifier=simplifier,
                                savefileprefix='',
                                copyfiletimes=copyfiletimes,
                                verbose=verbose,
                                verbosePeriod=verbosePeriod,
                                precision=precision,
                                restrictionList=restrictionList,
                                precomputedPartitionList=precomputedPartitionList,
                                solutionFixer=solutionFixer,
                                num_workers=num_workers,
                                jobqueuemaxsize=jobqueuemaxsize,
                                resultqueuemaxsize=resultqueuemaxsize)
        

#############################
#  Modernization functions  #
#############################

def modernizesolinsp(sol):
    """
    Takes a dictionary formatted (see solveFragBirkhoffDual) defrag solution sol whose 'insp' entry is in old
    format (i.e., either a set of strings or a dict with string keys, each string representing a partition) and
    updates it to the modern format. If sol does not have an 'insp' entry, then does nothing. If sol is already
    in modern format, then does nothing.
    """
    insp = sol.get('insp', None)
    if insp == None:
        return

    parameters = sol['parameters']

    lw = parameters.get('lightweight', None)

    if lw == None or (type(insp[-1]) is set) or lw:
        if lw == None or lw == True:
            parameters['lightweight'] = 1
        newinsp = []
        for entry in insp:
            if type(entry) is set: # set format
                newentry = set()
                for name in entry:
                    i = 0
                    while name[i] != '[':
                        i += 1
                    newentry.add((name[:i], int(name[i+1:-1])))
                newinsp.append(CompactSet.CompactSet(newentry))
            elif type(entry) is list: # compactify format
                newinsp.append(CompactSet.CompactSet((name, k)
                                                     for (r,a,b,s) in entry
                                                     for name in r
                                                     for k in range(a,b,s)))
            elif type(entry) is CompactSet.CompactSet:
                newinsp.append(entry)
            else:
                raise Exception('Format not recognized in entry:\n %s' % entry)
        sol['insp'] = tuple(newinsp)
    else:
        if lw == False:
            parameters['lightweight'] = 0
        newinsp = []
        for entry in insp:
            newentry = dict()
            for (name, val) in entry.items():
                if type(name) is tuple:
                    newentry = entry
                    break
                i = 0
                while name[i] != '[':
                    i += 1
                newentry[(name[:i], int(name[i+1:-1]))] = val
            newinsp.append(newentry)
        sol['insp'] = tuple(newinsp)
    

def modernizefile(filename, savefilename=None, *, copyfiletimes=True):
    """
    Reads a dictionary formatted solution (see solveFragBirkhoffDual) from file filename, applies
    modernizesolinsp to it and saves it in file savefilename.

    If savefilename is None (default), then saves on the same file from the input (overwritting it).

    If copyfiletimes is True, the access time and modified time of file filename will be copied to file 
    savefilename.
    """
    if savefilename == None:
        savefilename = filename

    sol = load(filename)
    if copyfiletimes:
        filestat = os.stat(filename)
    
    modernizesolinsp(sol)
    with open(savefilename, 'wb') as f:
        pickle.dump(sol, f, pickle.HIGHEST_PROTOCOL)

    if copyfiletimes:
        os.utime(savefilename, (filestat.st_atime, filestat.st_mtime))


def modernizeallfiles(foldername, extension='.pkl', *, copyfiletimes=True, verbose=True):
    """
    Applies modernizefile to all files from folder foldername (without trailing slash) with extension extension.
    The parameter savefilename from modernizefile used is the original file name.
    Parameter copyfiletimes is passed to modernizefile.
    If verbose is True, prints messages for monitoring.
    """
    for f in glob.glob(foldername + '/*' + extension):
        if verbose:
            print("Processing file %s. . ." % f)
        modernizefile(f, f, copyfiletimes=copyfiletimes)


##########################
#  Smart skip functions  #
##########################

def smartSkip(lmax, k0, *, skip=10):
    turnpoint = 2 * (lmax + 1)
    maxpoint = 2 * (lmax + k0)
    ret = list(range(2, turnpoint, 2))
    ret.extend(range(turnpoint+2, maxpoint+1, skip))
    return ret

def skipify(par, *, skip=10):
    for i in range(len(par)):
        (lmax, k0, c) = par[i]
        par[i] = (lmax, k0, c, smartSkip(lmax, k0, skip=skip))

def cleverSkip(lmax, k0, *, firstinc=4, smallskip=30, bigskip=50):
    turnpoint = 2 * (lmax + 1)
    maxpoint = 2 * (lmax + k0)
    ret = list(range(lmax+firstinc, turnpoint, smallskip))
    ret.extend(range(turnpoint+2, maxpoint+1, bigskip))
    return ret

def cskipify(par, *, firstinc=4, smallskip=30, bigskip=50):
    for i in range(len(par)):
        (lmax, k0, c) = par[i]
        par[i] = (lmax, k0, c, cleverSkip(lmax, k0, firstinc=firstinc, smallskip=smallskip, bigskip=bigskip))


##################################################
#              Deprecated functions              #
#  (_old_ is prepended to their original names)  #
##################################################

def _old_createBirkhoffPrimal(lmax, k0 = 15, c = fractions.Fraction(17,10), *, lset=None, mset=None,
                         kostkaRounder=dummy, negativeBelly=False, negativeLeg=False):
    """
    Creates the primal linear program associated to the Birkhoff problem.
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - kostkaRounder is a function to perform primal rounding of Kostka coefficients.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.

    Returns a tuple (obj, A, b, Aineq, vartype, obj0, maximization) to be used as parameters for
    simplex.solveNamedSimplex respectively as (c, A, b, Aineq, vartype, c0, maximization).
    """
    partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)]
    mmax = 2*k0 + 2*lmax

    A = []
    b = []
    obj = {}
    Aineq = []
    vartype = {}

    obj['Max'] = 1
    vartype['Max'] = 'u'

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))


    for l in lset:
        parsevaltable = Partition.parseval(l, k0, negativeBelly=negativeBelly, negativeLeg=negativeLeg)
        psiname = 'psi%03d' % l
        vartype[psiname] = 'u'

        A.append({psiname:1, 'Max':-1})
        Aineq.append('<=')
        b.append(0)

        rest = {}
        rest[psiname] = -1

        for bellysize in range(0, l+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    key = part.bname()
                    rest['%s[%03d]' % (key, k)] = 2 * parsevaltable[k].get(key, 0)

        A.append(rest)
        Aineq.append('=')
        # b.append(-2 + 2 * len(partitions[l]) * tailbound(l, k0, c))
        # The below uses joint bounds for the tails
        b.append(-2 + 2 * tailbound(l, k0, c))


    for m in mset:
        (Kostkakeys, Kostkatable) = Partition.kostka(m, lmax, k0)

        rest = {}

        for bellysize in range(0,lmax+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    key = part.bname()
                    rest['%s[%03d]' % (key, k)] = Kostkatable[k][key]

        rhs = kostkaRounder(rest, c ** m - 1)

        A.append(rest)
        Aineq.append('<=')
        b.append(rhs)

    for bellysize in range(0,lmax+1):
        for part in partitions[bellysize]:
            for k in range(max(part.height(),1),k0+1):
                vartype['%s[%03d]' % (part.bname(), k)] = 'p'

    return (obj, A, b, Aineq, vartype, 0, False)


def _old_createBirkhoffDual(lmax, k0 = 15, c = fractions.Fraction(17,10), *, lset=None, mset=None,
                       kostkaRounder=dummy, dualKostkaRounder=dualDummy, negativeBelly=False, negativeLeg=False):
    """
    Creates the dual linear program associated to the Birkhoff problem.
    - lmax is the maximum Parseval level used.
    - k0 is the truncation parameter.
    - c is the pseudorandomness parameter.
    - lset is the set of all Parseval levels used (if None, uses all even numbers from 0 to lmax, inclusive).
    - mset is the set of all Young restrictions used (if None, uses all numbers from 1 to 2*k0 + 2*lmax,
    inclusive).
    - kostkaRounder is a function to perform (dualization of the) primal rounding of Kostka coefficients.
    - dualKostkaRounder is a function to perform dual rounding of Kostka coefficients.
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.

    Returns a tuple (obj, A, b, Aineq, vartype, obj0, maximization) to be used as parameters for
    simplex.solveNamedSimplex respectively as (c, A, b, Aineq, vartype, c0, maximization).
    """

    partitions = [list(Partition.enumeratePartitions(l)) for l in range(lmax+1)]
    mmax = 2*k0 + 2*lmax

    A = []
    b = []
    obj = {}
    Aineq = []
    vartype = {}

    obj0 = 2

    if lset == None:
        lset = list(range(0,lmax+1,2))
    assert max(lset) <= lmax, 'lset (%s) has l > lmax = %s' % (str(lset), str(lmax))

    if mset == None:
        mset = list(range(1,mmax+1))
    assert max(mset) <= mmax, 'mset (%s) has m > mmax = %s' % (str(mset), str(mmax))

    for l in lset:
        # obj['w[%03d]' % l] = -2 * len(partitions[l]) * tailbound(l, k0, c)
        # The below uses joint bounds for the tails
        obj['w[%03d]' % l] = -2 * tailbound(l, k0, c)

        vartype['w[%03d]' % l] = 'p'
    for m in mset:
        vartype['y[%03d]' % m] = 'p'

    A.append(dict([('w[%03d]' % l, 1) for l in lset]))
    Aineq.append('=')
    b.append(1)

    restrictions = {(bellysize, part, k):{} for bellysize in range(0,lmax+1) for part in partitions[bellysize] for k in range(max(part.height(),1), k0+1)}
    
    for m in mset:
        Kostkatable = Partition.kostka(m,lmax,k0)[1]
        primalrest = {}
        for bellysize in range(0, lmax+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    primalrest[(k, part.bname())] = Kostkatable[k][part.bname()]

        primalrhs = kostkaRounder(primalrest, c ** m - 1)
        obj['y[%03d]' % m] = -primalrhs
        for bellysize in range(0, lmax+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    restrictions[(bellysize, part, k)]['y[%03d]' % m] = primalrest[(k,part.bname())]

    for l in lset:
        parsevaltable = Partition.parseval(l, k0, negativeBelly=negativeBelly, negativeLeg=negativeLeg)
        for bellysize in range(0, l+1):
            for part in partitions[bellysize]:
                for k in range(max(part.height(),1), k0+1):
                    restrictions[(bellysize, part, k)]['w[%03d]' % l] = 2 * parsevaltable[k].get(part.bname(), 0)

    for bellysize in range(0, lmax+1):
        for part in partitions[bellysize]:
            for k in range(max(part.height(),1), k0+1):
                A.append(restrictions[(bellysize, part, k)])
                Aineq.append('>=')
                b.append(0)

    obj0 = dualKostkaRounder(obj, A, b, Aineq, vartype, obj0, c=c)
                
    return (obj, A, b, Aineq, vartype, obj0, True)

def _old_modernizesolinsp(sol):
    """
    Takes a dictionary formatted (see solveFragBirkhoffDual) defrag solution sol whose 'insp' entry is in old
    format (i.e., either a set of strings or a dict with string keys, each string representing a partition) and
    updates it to the modern format. If sol does not have an 'insp' entry, then does nothing. If sol is already
    in modern format, then does nothing.
    """
    insp = sol.get('insp', None)
    if insp == None:
        return

    parameters = sol['parameters']

    if parameters.get('lightweight', False) or type(insp[-1]) is set:
        parameters['lightweight'] = True
        newinsp = []
        for entry in insp:
            if type(entry) is set:
                newentry = set()
                for name in entry:
                    i = 0
                    while name[i] != '[':
                        i += 1
                    newentry.add((name[:i], int(name[i+1:-1])))
                newinsp.append(compactify(newentry))
            else:
                newinsp.append(entry)
        sol['insp'] = tuple(newinsp)
    else:
        parameters['lightweight'] = False
        newinsp = []
        for entry in insp:
            newentry = dict()
            for (name, val) in entry.items():
                if type(name) is tuple:
                    newentry = entry
                    break
                i = 0
                while name[i] != '[':
                    i += 1
                newentry[(name[:i], int(name[i+1:-1]))] = val
            newinsp.append(newentry)
        sol['insp'] = tuple(newinsp)

#################################################################
#  Compact sets (deprecated, we now use the module CompactSet)  #
#################################################################

def _old_compactify(s):
    """
    Given a set s of pairs (name, k) where name is a string and k is an int, returns a compact set
    representing it (some reasonable attempt at a short representation is done, no optimality or near
    optimality guarantees are given).

    A compact set is a list whose entries are quadruples (r, a, b, s), where r is a set of strings and a, b, s
    are ints. The quadruple (r, a, b, s) represents the set set((name, k) for name in r for k in range(a,b,s))
    and the compact set represents the union of the sets represented by its elements. See also decompactify.
    """
    d = dict()
    for (name, k) in s:
        d.setdefault(name, []).append(k)

    d2 = dict()
    for (name, klist) in d.items():
        klist.sort()
        start = klist[0]
        step = None
        for i in range(1, len(klist)):
            if step == None:
                step = klist[i] - klist[i-1]
            elif klist[i] != klist[i-1] + step:
                d2.setdefault((start, klist[i-1] + 1, step), set()).add(name)
                start = klist[i]
                step = None
        if step == None:
            step = 1
        d2.setdefault((start, klist[-1] + 1, step), set()).add(name)
            
    return [(r, a, b, s) for ((a,b,s),r) in d2.items()]

def _old_decompactify(s):
    """
    Given a compact set s, returns the set it represents.

    See compactify for the definition of compact sets.
    """
    return set((name, k) for (r,a,b,s) in s for name in r for k in range(a,b,s))
