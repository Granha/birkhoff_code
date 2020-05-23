def withoutConsecutiveRepetitions(l):
    """
    Returns a generator that iterates over l removing consecutive repetitions (keeps only first copy).
    """
    it = iter(l)
    try:
        elem = next(it)
    except StopIteration:
        # l is empty
        return
    try:
        while True:
            newelem = next(it)
            if newelem != elem:
                yield elem
                elem = newelem
    except StopIteration:
        yield elem

class CompactIntegerSet:
    """
    A class supporting element addition for a reasonably compact representation of sets of int.

    The entry cis of a CompactIntegerSet is a list of triples (a, b, s) of ints.
    A triple (a, b, s) represents the set set(range(a,b,s)) and the CompactIntegerSet represents the
    union of the sets represented by all triples in cis.

    There is no guarantee that each int is uniquely represented in the union above.

    Some reasonable attempt at a short representation is done, no optimality or near optimality guarantees are
    given.

    The following invariants are mantained.
    - All triples (a, b, s) have s >= 1.
    - For a triple (a, b, s) the last element of range(a, b, s) is b - 1 (even if s is not 1).
    - If a triple (a, b, s) is such that range(a, b, s) has only one element, then s = 1.

    See also SemiCompactSet and CompactSet.
    """

    def __init__(self, other=[], *, sortBefore=True, reverseTraversal=True, lazyLookup=True):
        """
        Initializes self from either a CompactIntegerSet or an iterable that provides pairs integers.

        In the case of an iterable that is not a CompactIntegerSet, uses the function update
        and passes on the parameters to sortBefore, reverseTraversal and lazyLookup. To get this behaviour
        for CompactIntegerSet use CompactIntegerSet((x for x in other)).
        """
        if isinstance(other, CompactIntegerSet):
            self.cis = other.cis.copy()
        else:
            self.cis = list()
            self.update(other, sortBefore=sortBefore, reverseTraversal=reverseTraversal, lazyLookup=lazyLookup)

    def add(self, k, *, reverseTraversal=True, lazyLookup=False):
        """
        Adds int k to self.

        If reverseTraversal is True, looks for appropriate place for insertion in reverse order (this makes
        adding elements faster if several additions are made in monotonic order).

        If lazyLookup is True, makes only a lazy attempt at compactification by looking up only one triple
        (this makes addition faster but may not compactify the set).

        See also update.
        """
        length = len(self.cis)
        if length > 0 and lazyLookup:
            length = 1
        (start, pastend, step) = (-1, -length-1, -1) if reverseTraversal else (0, length, 1)

        for i in range(start, pastend, step):
            (a, b, s) = self.cis[i]

            if b - 1 + s == k:
                # k allows us to extend the triple by one more integer forward
                self.cis[i] = (a, k + 1, s)
                return

            if a - s == k:
                # k allows us to extend the triple by one more integer backward
                self.cis[i] = (k, b, s)
                return

            if b - 1 == a:
                # The triple represents only one element
                if a < k:
                    # k can be added to it forward
                    self.cis[i] = (a, k + 1, k - a)
                elif k < a:
                    # k can be added to it backward
                    self.cis[i] = (k, a+1, a - k)
                # Otherwise, this unique element is k
                
                return

            if a <= k and k < b and (k - a) % s == 0:
                # (name, k) is already in self
                return

        # k cannot be added by extending any triple
        self.cis.append((k, k+1, 1))

        # Warning: it may seem that the code above avoids an element from being represented more than once,
        # but this is not the case: consider the following addition order (with reverseTraversal=True).
        # 1 -> [(1,2,1)]
        # 2 -> [(1,3,1)]
        # 4 -> [(1,3,1), (4,5,1)]
        # 3 -> [(1,3,1), (3,5,1)]
        # 2 -> [(1,3,1), (2,5,1)] (2 is represented twice)

    def update(self, elems, *, sortBefore=True, reverseTraversal=True, lazyLookup=False):
        """
        Adds all elements of an iterable elems of int to self.

        The final result represents the same set as if add was called for each element of elems.

        The parameters reverseTraversal and lazyLookup are passed on to add.

        If sortBefore is True, sorts the integers and removes repetitions before adding them (this along
        with reverseTraversal=True and lazyLookup=True can potentially make this call faster than calling
        add for each element of elems).
        """

        if sortBefore:
            elems = sorted(elems)

        for e in elems:
            self.add(e, reverseTraversal=reverseTraversal, lazyLookup=lazyLookup)

    def __iter__(self):
        """
        Returns an iterator to the elements in the set represented by self (possibly with repetitions).
        """
        for (a,b,s) in self.cis:
            for k in range(a,b,s):
                yield k

    def copy(self):
        """
        Returns a copy of self.
        """
        return CompactIntegerSet(self)

    def clear(self):
        """
        Removes all elements of self.
        """
        self.cis = list()

    def reorganize(self):
        """
        Removes all repetitions and makes another attempt at reducing the representation size.

        Warning: this function traverses the whole container, so it may be slow.
        """
        l = list(self)
        self.clear()

        self.update(l, sortBefore=True, reverseTraversal=True, lazyLookup=True)

    def __repr__(self):
        """
        String representation.
        """
        return self.cis.__repr__()

    def __len__(self):
        """
        Length of representation (this is not the size of the set).
        It is 0 if and only if self represents the empty set.
        """
        return len(self.cis)


class SemiCompactSet:
    """
    A class supporting element addition for a reasonably compact representation of sets of pairs (name, k),
    where name is hashable and k is an int.

    The entry scs of a SemiCompactSet is a dict whose entries are name:r, where name is hashable and r is a
    CompactIntegerSet. An entry name:r represents the set
    set((name, k) for k in r)
    and the SemiCompactSet represents the union of the sets represented by all elements of scs.

    There is no guarantee that each pair (name, k) is uniquely represented in the union above.

    Some reasonable attempt at a short representation is done, no optimality or near optimality guarantees are
    given.

    See also CompactSet.
    """
    
    def __init__(self, other=[], *, sortBefore=True, reverseTraversal=True, lazyLookup=True):
        """
        Initializes self from either a SemiCompactSet, a CompactSet or an iterable that provides pairs 
        (name, k), where name is hashable and k is an int.

        In the case of a CompactSet or an iterable that is not a SemiCompactSet, uses the function update
        and passes on the parameters to sortBefore, reverseTraversal and lazyLookup. To get this behaviour
        for SemiCompactSet or CompactSet use SemiCompactSet((x for x in other)).
        """
        if isinstance(other, SemiCompactSet):
            self.scs = {name:r.copy() for (name, r) in other.scs.items()}
        elif isinstance(other, CompactSet):
            self.scs = dict()
            for (r,a,b,s) in other.cs:
                for name in r:
                    self.scs.setdefault(name, CompactIntegerSet()).update(range(a,b,s),
                                                                          sortBefore=sortBefore,
                                                                          reverseTraversal=reverseTraversal,
                                                                          lazyLookup=lazyLookup)
        else:
            self.scs = dict()
            self.update(other, sortBefore=sortBefore, reverseTraversal=reverseTraversal, lazyLookup=lazyLookup)

    def add(self, elem, *, reverseTraversal=True, lazyLookup=False):
        """
        Adds element elem (assumed to be a pair (name,k), where name is hashable and k is an int) to self.

        Passes reverseTraversal and lazyLookup to CompactIntegerSet.add.

        See also update.
        """
        (name, k) = elem
        self.scs.setdefault(name, CompactIntegerSet()).add(k,
                                                           reverseTraversal=reverseTraversal,
                                                           lazyLookup=lazyLookup)

    def update(self, elems, *, sortBefore=True, reverseTraversal=True, lazyLookup=False):
        """
        Adds all elements of an iterable elems (all elements are assumed to be pairs (name,k), where name is
        hashable and k is an int) to self.

        The final result represents the same set as if add was called for each element of elems.

        The parameters sortBefore, reverseTraversal and lazyLookup are passed on to CompactIntegerSet.update.
        """

        d = dict()

        for (name, k) in elems:
            d.setdefault(name, list()).append(k)

        for (name, klist) in d.items():
            self.scs.setdefault(name, CompactIntegerSet()).update(klist,
                                                                  sortBefore=sortBefore,
                                                                  reverseTraversal=reverseTraversal,
                                                                  lazyLookup=lazyLookup)

    def __iter__(self):
        """
        Returns an iterator to the elements in the set represented by self (possibly with repetitions).
        """
        for (name, r) in self.scs.items():
            for k in r:
                yield (name, k)

    def copy(self):
        """
        Returns a copy of self.
        """
        return SemiCompactSet(self)

    def clear(self):
        """
        Removes all elements of self.
        """
        self.scs = dict()

    def reorganize(self):
        """
        Removes all repetitions and makes another attempt at reducing the representation size.

        Warning: this function traverses the whole container, so it may be slow.
        """
        for (name, r) in self.scs.items():
            r.reorganize()

    def __repr__(self):
        """
        String representation.
        """
        return self.scs.__repr__()

    def __len__(self):
        """
        Length of dictionary of the representation (this is not the size of the set).
        It is 0 if and only if self represents the empty set.
        """
        return len(self.scs)


class CompactSet:
    """
    A class for a reasonably compact representation of sets of pairs (name, k), where name is hashable and
    k is an int.

    The entry cs of a CompactSet is a list whose entries are quadruples (r, a, b, s), where r is a set of
    hashables. The quadruple (r, a, b, s) represents the set
    set((name, k) for name in r for k in range(a,b,s))
    and the CompactSet represents the union of the sets represented by its elements.

    There is no guarantee that each pair (name, k) is uniquely represented in the union above.

    Some reasonable attempt at a short representation is done, no optimality or near optimality guarantees are
    given.

    The following invariants are mantained.
    - All quadruples (r, a, b, s) have s >= 1.
    - For a quadruple (r, a, b, s) the last element of range(a, b, s) is b - 1 (even if s is not 1).
    - If a quadruple (r, a, b, s) is such that range(a, b, s) has only one element, then s = 1.

    See also SemiCompactSet.
    """

    def __init__(self, other=[]):
        """
        Initializes self from either a CompactSet, a SemiCompactSet or an iterable that provides pairs
        (name, k), where name is hashable and k is an int.
        """
        if isinstance(other, CompactSet):
            self.cs = [(r.copy(), a, b, s) for (r, a, b, s) in other.cs]
            return

        if not isinstance(other, SemiCompactSet):
            other = SemiCompactSet(other)

        d = dict()
        for (name, cis) in other.scs.items():
            for krange in cis.cis:
                d.setdefault(krange, set()).add(name)
        self.cs = [(r, a, b, s) for ((a, b, s), r) in d.items()]

    def __iter__(self):
        """
        Returns an iterator to the elements in the set represented by self (possibly with repetitions).
        """
        for (r, a, b, s) in self.cs:
            for name in r:
                for k in range(a, b, s):
                    yield (name, k)

    def copy():
        """
        Returns a copy of self.
        """
        return CompactSet(self)

    def clear(self):
        """
        Removes all elements of self.
        """
        self.cs = list()

    def reorganize(self):
        """
        Removes all repetitions and makes another attempt at reducing the representation size.

        Warning: this function traverses the whole container, so it may be slow.
        """
        scs = SemiCompactSet(self)
        scs.reorganize()
        self.__init__(scs)

    def __repr__(self):
        """
        String representation.
        """
        return self.cs.__repr__()

    def __len__(self):
        """
        Length of the representation (this is not the size of the set).
        It is 0 if and only if self represents the empty set.
        """
        return len(self.cs)
