import unittest

from Partition import *

class TestPartition(unittest.TestCase):

    def test_height(self):
        N = 100
        layers = []
        for i in range(1,N+1):
            layers.append(N)
            N -= 1
            p = Partition(layers)
            self.assertEqual(p.height(), i)

    def test_basic(self):
        p = Partition([5,5,3,3,1])
        self.assertEqual(p.size(), 17)
        self.assertEqual(p.empty(), False)
        self.assertEqual(p.height(), 5)
        self.assertEqual(p.width(), 5)
        self.assertEqual(p.transpose().layers, [5,4,4,2,2])
        self.assertEqual(p.findBelly().layers, [4,2,2])
        self.assertEqual(p.underlinePrint(), '5_5_3_3_1')
        self.assertEqual(p.bname(), 'b5_5_3_3_1')

        jp = JointPartition(17)
        self.assertEqual(jp.size(), 17)
        self.assertEqual(jp.transpose(), jp)
        self.assertEqual(jp.underlinePrint(), '{17}')
        self.assertEqual(jp.bname(), 'b{17}')
        
    def test_dimension(self):
        p = Partition([5,5,3,3,1])
        self.assertEqual(p.dimension(), 1361360)

    def test_dimension2(self):
        partitions = [Partition([]),
                      Partition([1]),
                      Partition([2]), Partition([1,1]),
                      Partition([3]), Partition([2,1]), Partition([1,1,1]),
                      Partition([4]), Partition([3,1]), Partition([2,2]), Partition([2,1,1]), Partition([1,1,1,1]),
                      #
                      Partition([5]), Partition([4,1]), Partition([3,2]), Partition([3,1,1]), Partition([2,2,1]),
                      Partition([2,1,1,1]), Partition([1,1,1,1,1]),
                      #
                      Partition([6]), Partition([5,1]), Partition([4,2]), Partition([4,1,1]), Partition([3,3]),
                      Partition([3,2,1]), Partition([3,1,1,1]), Partition([2,2,2]), Partition([2,2,1,1]),
                      Partition([2,1,1,1,1]), Partition([1,1,1,1,1,1])]
        dimensions = [1,
                      1,
                      1,1,
                      1,2,1,
                      1,3,2,3,1,
                      1,4,5,6,5,4,1,
                      1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]
        for i in range(len(partitions)):
            self.assertEqual(partitions[i].dimension(), dimensions[i])

    def test_identities(self):
        p = Partition([7,5,3,3,1])
        self.assertEqual(p.height(), p.transpose().width())
        self.assertEqual(p.width(), p.transpose().height())
        self.assertEqual(p.transpose().transpose(), p)
        self.assertEqual(p.transpose().findBelly(), p.findBelly().transpose())
        self.assertEqual(p.size(), p.height() + p.width() + p.findBelly().size() - 1)
        self.assertEqual(bellyPartition(p.size(), p.height() - 1, p.findBelly()), p)
        self.assertEqual(addColumn(p.height() - 1, p.findBelly()).layers, p.layers[1:])
        self.assertEqual(p.dimension(), p.transpose().dimension())

    def test_enumeratePartitions(self):
        correct = [   1,    1,    2,    3,    5,    7,   11,   15,   22,   30,
                     42,   56,   77,  101,  135,  176,  231,  297,  385,  490,
                    627,  792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565,
                   5604] # from 0 to 30

        for l in range(len(correct)):
            self.assertEqual(len(list(enumeratePartitions(l))), correct[l])

    def test_enumerateArmParseval(self):
        partitions6 = [Partition([6]), Partition([5,1]), Partition([4,2]), Partition([4,1,1]), Partition([3,3]),
                       Partition([3,2,1]), Partition([3,1,1,1]), Partition([2,2,2]), Partition([2,2,1,1]),
                       Partition([2,1,1,1,1]), Partition([1,1,1,1,1,1])]
        dimensions6 = [1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]
        heights6 = [1, 2, 2, 3, 2, 3, 4, 3, 4, 5, 6]
        bellies6 = [Partition([]), Partition([]), Partition([1]), Partition([]), Partition([2]),
                    Partition([1]), Partition([]), Partition([1,1]), Partition([1]),
                    Partition([]), Partition([])]

        correct = [(heights6[i] - 1, bellies6[i], dimensions6[i]) for i in range(len(partitions6)) if heights6[i] >= 2]
        
        self.assertEqual(list(enumerateArmParseval(6)), correct)

    def test_bellyOfAddBorder(self):
        p = Partition([5,5,3,3,1])
        added = [Partition([6,5,3,3,1]), Partition([6,6,3,3,1]), Partition([6,6,6,3,1]), Partition([6,6,6,4,1]),
                 Partition([6,6,6,4,4])]
        bellyOfAdded = [Partition([4,2,2]), Partition([5,2,2]), Partition([5,5,2]), Partition([5,5,3]),
                        Partition([5,5,3,3])]
        for i in range(p.height()):
            result = bellyOfAddBorder(p,i)
            self.assertEqual(result, added[i].findBelly())
            self.assertEqual(result, bellyOfAdded[i])

    def test_enumerateBellyParseval(self):
        correct = [(1, Partition([5]), -5), # dim of [5,1]
                   (1, Partition([4]), -9), # dim of [4,2]
                   (2, Partition([4]), -10), (2, Partition([4,1]), +10), # dim of [4,1,1]
                   (1, Partition([3]), -5), # dim of [3,3]
                   (2, Partition([3]), -16), (2, Partition([3,2]), +16), # dim of [3,2,1]
                   (3, Partition([3]), -10), (3, Partition([3,1]), +10), (3, Partition([3,1,1]), -10),
                   # dim of [3,1,1,1]
                   (2, Partition([2,1]), -5), (2, Partition([2,2]), +5), # dim of [2,2,2]
                   (3, Partition([2]), -9), (3, Partition([2,2]), +9), (3, Partition([2,2,1]), -9),
                   # dim of [2,2,1,1]
                   (4, Partition([2]), -5), (4, Partition([2,1]), +5), (4, Partition([2,1,1]), -5),
                   (4, Partition([2,1,1,1]), +5),
                   # dim of [2,1,1,1,1]
                   (5, Partition([1]), -1), (5, Partition([1,1]), +1), (5, Partition([1,1,1]), -1),
                   (5, Partition([1,1,1,1]), +1), (5, Partition([1,1,1,1,1]), -1)
                   # dim of [1,1,1,1,1,1]
                   ]

        self.assertEqual(set(enumerateBellyParseval(6)), set(correct))

    def test_enumerateLegParseval(self):
        partitions6 = [Partition([6]), Partition([5,1]), Partition([4,2]), Partition([4,1,1]), Partition([3,3]),
                       Partition([3,2,1]), Partition([3,1,1,1]), Partition([2,2,2]), Partition([2,2,1,1]),
                       Partition([2,1,1,1,1]), Partition([1,1,1,1,1,1])]
        dimensions6 = [1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1]
        heights6 = [1, 2, 2, 3, 2, 3, 4, 3, 4, 5, 6]
        k0 = 11

        correct = [(k, partitions6[i], (-1)**k * dimensions6[i]) for i in range(len(partitions6)) for k in range(max(heights6[i],1),k0+1)]

        self.assertEqual(set(enumerateLegParseval(l=6, k0=k0)), set(correct))

    def test_parseval_singlePartitionParseval(self):
        l = 6
        k0 = 11

        partitions = [list(enumeratePartitions(i)) for i in range(l+1)]

        for nB in [False, True]:
            for nL in [False, True]:
                pars = parseval(l, k0, negativeBelly=nB, negativeLeg=nL, precomputedPartitionList=partitions[l])

                for i in range(l+1):
                    for belly in partitions[i]:
                        for k in range(max(belly.height(), 1), k0+1):
                            self.assertEqual(pars[k].get(belly.bname(), 0),
                                             singlePartitionParseval(l, k, belly,
                                                                     negativeBelly=nB, negativeLeg=nL),
                                             '\nFailed with %s[%03d]\n' % (belly.bname(), k)
                                             + 'negativeBelly=%d\n' % nB
                                             + 'negativeLeg=%d\n' % nL)

    def test_parseval_singlePartitionParseval2(self):
        l = 10
        k0 = 20

        partitions = [list(enumeratePartitions(i)) for i in range(l+1)]

        for nB in [False, True]:
            for nL in [False, True]:
                pars = parseval(l, k0, negativeBelly=nB, negativeLeg=nL, precomputedPartitionList=partitions[l])

                for i in range(l+1):
                    for belly in partitions[i]:
                        for k in range(max(belly.height(), 1), k0+1):
                            self.assertEqual(pars[k].get(belly.bname(), 0),
                                             singlePartitionParseval(l, k, belly,
                                                                     negativeBelly=nB, negativeLeg=nL),
                                             '\nFailed with %s[%03d]\n' % (belly.bname(), k)
                                             + 'negativeBelly=%d\n' % nB
                                             + 'negativeLeg=%d\n' % nL)                            

if __name__ == '__main__':
    unittest.main()
