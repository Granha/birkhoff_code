from math import factorial
from functools import reduce

###########################################
#  Partition class and related functions  #
###########################################

class Partition:
    """
    Class for partitions.
    """

    def __init__(self, t):
        """
        Initializes partition with iterable t.
        
        Assumption: t represents a partition.
        """
        self.layers = list(t)

    def __repr__(self):
        """
        String representation.
        """
        return self.layers.__repr__()

    def size(self):
        """
        Returns size of the partition.
        """
        return sum(self.layers)

    def empty(self):
        """
        Returns True if empty partition.
        """
        return len(self.layers) == 0

    def height(self):
        """
        Returns height of the partition.
        """
        return len(self.layers)

    def width(self):
        """
        Returns width of the partition.
        """
        return self.layers[0] if not self.empty() else 0

    def transpose(self):
        """
        Returns transpose of the partition.
        """

        r = [0] * self.width()
        for v in self.layers:
            for i in range(v):
                r[i] += 1
        return Partition(r)

    def dimension(self):
        """
        Returns the dimension of the Specht Module associated with the partition.
        Computed via hook-length formula.
        """
        t = self.transpose()
        r = factorial(self.size())
        for i in range(len(self.layers)):
            for j in range(self.layers[i]):
                r //= (self.layers[i] - j) + (t.layers[j] - i) - 1
        return r

    def findBelly(self):
        """
        Returns the belly of the partition.
        """
        return Partition([v-1 for v in self.layers[1:] if v >= 2])

    def underlinePrint(self):
        """
        Returns a representation of the partition with underlines separating layers.
        """
        return '_'.join((str(v) for v in self.layers))

    def bname(self):
        """
        Returns the bname of the partition.
        """
        return 'b' + self.underlinePrint()

    def __eq__(self, other):
        """
        Returns True if partition is equal to other and False if they are different.
        (other is assumed to be of type Partition.)
        """
        return self.layers == other.layers

    def __hash__(self):
        """
        Returns hash value.
        """
        return hash(repr(self))

def frombname(s):
    """
    Returns the Partition represented by the string s, assumed to be in bname format.
    """
    
    return Partition([int(sl) for sl in s[1:].split('_') if sl != ''])

def bellyPartition(n, k, p):
    """
    Returns partition $b_{k,p}^n$ given by belly notation.
    """
    return Partition([n - k - p.size()] + [v+1 for v in p.layers] + [1] * (k - p.height()))

def addColumn(k, p):
    """
    Returns the partition obtained from p by adding a column of height k to the left.

    Assumption k >= p.height().
    """
    return Partition([v+1 for v in p.layers] + [1] * (k - p.height()))

def enumeratePartitions(n):
    """
    Returns a generator of all partitions of size n.
    """
    if (n == 0):
        yield Partition([])
        return
    r = (n,)
    yield Partition(r)
    while True:
        i = len(r) - 1
        while i > -1:
            if r[i] > 1:
                break
            i -= 1
        if i == -1:
            return
        s = len(r)-i
        r = r[:i] + (r[i] - 1,)
        while s > 0:
            r += (min(r[-1], s),)
            s -= r[-1]
        yield Partition(r)

def countPartitions(n, *, width=None):
    """
    Returns the number of partitions of n.
    If width is not None, returns the partitions of n of width width instead.
    """
    if width == None:
        return sum((countPartitions(n, width=i) for i in range(n+1)))
    if width > n:
        return 0
    if n == 0:
        return 1
    if width == 0:
        return 0

    tblen = len(countPartitions.table)
    while tblen < n:
        tblen += 1
        countPartitions.table.append([0] * tblen)

    if countPartitions.table[n-1][width-1] == 0:
        countPartitions.table[n-1][width-1] = sum((countPartitions(n-width, width=i) for i in range(width+1)))

    return countPartitions.table[n-1][width-1]

countPartitions.table = []



################################################
#  JointPartition class and related functions  #
################################################

class JointPartition:
    """
    A class for representing all partitions of a certain size.

    All functions that would be well-defined of the class Partition.Partition are implemented here,
    e.g., the transpose of an instance of JointPartition is (a copy of) itself, but the function height is
    not well-defined, so it is not implemented. Furthermore, string representation and related functions are
    also defined. This is done to allow usage of a JointPartition instance in place of a Partition instance
    in functions that only call members that are well-defined for JointPartition.
    """

    def __init__(self, s):
        """
        Initializes JointPartition to represent size s.
        """
        self.s = s

    def __repr__(self):
        """
        String representation.
        """
        return self.bname()

    def size(self):
        """
        Returns size of partitions represented by self.
        """
        return self.s

    def empty(self):
        """
        Returns True if class represents the class of the empty partition.
        """
        return self.s == 0

    def transpose(self):
        """
        Returns the JointPartition instance representing transpose of self (which is a copy of self).
        """
        return JointPartition(self.s)

    def underlinePrint(self):
        """
        Returns a representation of the partition class.
        """
        return '{%d}' % self.s

    def bname(self):
        """
        Returns the bname of the partition class.
        """
        return 'b' + self.underlinePrint()

    def __eq__(self, other):
        """
        Returns True if partition class is equal to other and False if they are different.
        (other is assumed to be of type JointPartition.)
        """
        return self.s == other.s

    def __hash__(self):
        """
        Returns hash value.
        """
        return hash(repr(self))


################################
#  Parseval related functions  #
################################

def enumerateArmParseval(l, *, precomputedPartitionList=None):
    """
    Returns a generator of triples (k,b,m) where (k,b) are the partitions in belly notation after removal
    of n-l completely contained in the arm and m is the multiplicity (with sign).

    - precomputedPartitionList is assumed to be one of the following.
    -- A list of all partitions of size exactly l.
    -- None, in which case such list will be computed by this function.
    """
    lpartitions = enumeratePartitions(l) if precomputedPartitionList == None else precomputedPartitionList
    
    return ((p.height()-1, p.findBelly(), p.dimension()) for p in lpartitions if p.height() >= 2)

def bellyOfAddBorder(p, i):
    """
    Adds a border that ends at line i to the partition p and returns the belly of the resulting partition.

    Assumption: 0 <= i <= p.height()-1.
    """
    return Partition([v for v in p.layers[:i]] + [v-1 for v in p.layers[i+1:] if v >= 2])

def enumerateBellyParseval(l, *, precomputedPartitionList=None):
    """
    Returns a generator of triples (k,b,m) where (k,b) are the partitions in belly notation after removal
    of n-l entering the belly (but not leaving it) and m is the multiplicity (with sign).

    - precomputedPartitionList is assumed to be one of the following.
    -- A list of all partitions of size exactly l.
    -- None, in which case such list will be computed by this function.
    """
    lpartitions = enumeratePartitions(l) if precomputedPartitionList == None else precomputedPartitionList

    return ((p.height()-1, bellyOfAddBorder(p,i), (-1)**i * p.dimension()) for p in lpartitions for i in range(1, p.height()))

def enumerateLegParseval(l, k0, *, precomputedPartitionList=None):
    """
    Returns a generator of triples (k,b,m) where (k,b) are the partitions in belly notation after removal
    of n-l containing both the arm and the leg and m is the multiplicity (with sign).

    - precomputedPartitionList is assumed to be one of the following.
    -- A list of all partitions of size exactly l.
    -- None, in which case such list will be computed by this function.
    """
    lpartitions = enumeratePartitions(l) if precomputedPartitionList == None else precomputedPartitionList
    
    return ((k, b, (-1)**k * b.dimension()) for b in lpartitions for k in range(max(b.height(),1),k0+1))

def parseval(l, k0, *, negativeBelly=False, negativeLeg=False, precomputedPartitionList=None):
    """
    Returns a list ret of Parseval coefficients corresponding to psi_l with truncation parameter k0.
    ret[k][name] is the Parseval coefficient of partition with belly with bname name and leg k.
    Omitted entries are assumed to be zero (but there may be explicitly zero entries).
    (Returns an empty list if l is odd.)
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.
    - precomputedPartitionList is assumed to be one of the following.
    -- A list of all partitions of size exactly l.
    -- None, in which case such list will be computed by this function.

    Assumes k0 >= l.
    """
    lpartitions = list(enumeratePartitions(l)) if precomputedPartitionList == None else precomputedPartitionList
    
    ret = [dict() for k in range(k0+1)]

    if l % 2 == 1:
        return []

    for (k, b, m) in enumerateArmParseval(l, precomputedPartitionList=lpartitions):
        ret[k][b.bname()] = m

    for (k, b, m) in enumerateBellyParseval(l, precomputedPartitionList=lpartitions):
        if negativeBelly and m > 0:
            m = -m
        ret[k][b.bname()] = m

    for (k, b, m) in enumerateLegParseval(l, k0, precomputedPartitionList=lpartitions):
        if negativeLeg and m > 0:
            m = -m
        ret[k][b.bname()] = m

    return ret

def singlePartitionParseval(l, k, belly, *, negativeBelly=False, negativeLeg=False):
    """
    Returns the Parseval coefficient corresponding to psi_l and the partition with leg k and belly belly.
    (Returns 0 if l is odd.)
    - if negativeBelly is True, changes all Parsevals of shapes with removals ending inside bellies to be negative.
    - if negativeLeg is True, changes all Parsevals of shapes with removals ending on legs to be negative.

    Assumption: k >= 1.
    """

    if l % 2 == 1:
        return 0

    if belly.size() > l:
        # No removal possible
        return 0

    if belly.size() == l:
        # Removal takes full border
        sgn = -1 if negativeLeg else (-1) ** k
        return sgn * belly.dimension() 

    remainPart = Partition([l + belly.width()] + [v+1 for v in belly.layers] + [1] * (k - belly.height()))
    remainSize = remainPart.size()

    if remainSize - remainPart.layers[0] + remainPart.layers[1] <= l:
        # Removal ends on arm
        remainPart.layers[0] -= remainSize - l
        return remainPart.dimension()


    # Appending auxiliary 0 at the end
    remainPart.layers.append(0)

    for i in range(belly.height()):
        remainSize -= remainPart.layers[i] - remainPart.layers[i+1] + 1
        if remainSize < l:
            # No removal is possible at all
            return 0
        remainPart.layers[i] = remainPart.layers[i+1] - 1

        if remainPart.layers[i+2] == remainPart.layers[i+1]: # Auxiliary 0 is needed for this
            # It is not possible to stop at this level as next level has same size
            continue

        if remainSize - remainPart.layers[i+1] + remainPart.layers[i+2] <= l:
            # Removal ends at this level
            if remainSize == l:
                # Removal would require no cell to be removed at this level but leave previous level with 
                # size equal to this level minus 1.
                return 0
            
            remainPart.layers[i+1] -= remainSize - l
            remainPart.layers = remainPart.layers[:-1] # Removing auxiliary 0
            sgn = -1 if negativeBelly else (-1) ** (i+1)
            return sgn * remainPart.dimension()

    # If this point is reached, then no removal is possible at all (it would end strictly inside the leg)
    return 0


##############################
#  Kostka related functions  #
##############################

def fallingFactorial(n,k):
    """
    Returns $(n)_k$.
    """
    return reduce(lambda i,j : i*j, (n-i for i in range(k)), 1)

def binomial(n,k):
    """
    Returns $\binom{n}{k}$.
    """
    return fallingFactorial(n,k) // factorial(k)

def kostkaBinomial(m, k, bellysize):
    """
    Returns the binomial of the Kostka number for partition with belly size bellysize and leg k corresponding
    to the Young module of the hook m.
    """
    return binomial(m, k + bellysize)

def kostkaDimension(k, belly):
    """
    Returns the Kostka dimension for partition with belly belly and leg k corresponding to the Young module
    of the hooks.
    """
    return addColumn(k, belly).dimension()

def kostkaNumber(m, k, belly):
    """
    Returns the Kostka number for partition with belly belly and leg k corresponding to the Young module
    of the hook m.
    """
    return kostkaBinomial(m, k, belly.size()) * kostkaDimension(k, belly)


##################################################
#              Deprecated functions              #
#  (_old_ is prepended to their original names)  #
##################################################

def _old_parseval(l, k0, *, negativeBelly=False, negativeLeg=False):
    """
    Return a pair (names, a), where names is the list of partition bnames and a is the array of
    Parseval coefficients corresponding to psi_l. (Returns a zero table if l is odd.)

    Assumes k0 >= l.
    """ 
    partitions = sum((list(enumeratePartitions(i)) for i in range(l+1)), [])
    printpart = [p.bname() for p in partitions]
    ret = [dict(zip(printpart, [0]*len(printpart))) for k in range(k0+1)]

    if l % 2 == 1:
        return (printpart, ret)

    for (k, b, m) in enumerateArmParseval(l):
        ret[k][b.bname()] = m

    for (k, b, m) in enumerateBellyParseval(l):
        if negativeBelly and m > 0:
            m = -m
        ret[k][b.bname()] = m

    for (k, b, m) in enumerateLegParseval(l, k0):
        if negativeLeg and m > 0:
            m = -m
        ret[k][b.bname()] = m

    return (printpart, ret)

def kostka(m, lmax, k0):
    """
    Return a pair (names, a), where names is the list of partition bnames and a is the array of
    Kostka numbers corresponding to the Young module of the hook m.
    """
    partitions = sum((list(enumeratePartitions(i)) for i in range(lmax+1)), [])
    printpart = [p.bname() for p in partitions]
    ret = [dict(zip(printpart, [0]*len(printpart))) for k in range(k0+1)]

    for k in range(1,k0+1):
        for b in partitions:
            if k >= b.height():
                ret[k][b.bname()] = binomial(m, k + b.size()) * addColumn(k, b).dimension()

    return (printpart, ret)



# Local Variables:
# mode: python
# mode: visual-line
# End:
