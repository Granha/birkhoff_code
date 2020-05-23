import unittest

from simplex import *
from fractions import Fraction

class TestSimplex(unittest.TestCase):

    def test_program(self):
        c = [1,2,3,4]
        A = [[1,-1,3,1],[-1,2,-1,-1],[2,1,-7,1]]
        b = [5,-2,3]

        c = [Fraction(x) for x in c]
        A = [[Fraction(x) for x in line] for line in A]
        b = [Fraction(x) for x in b]

        solution = [0, Fraction(13, 7), Fraction(4, 7), Fraction(36, 7)]
        value = sum(solution[var] * c[var] for var in range(len(c)))

        self.assertEqual(solveStandardSimplex(c=c,A=A,b=b),(solution,value))

    def test_program2(self):
        c = [1,1]
        A = [[1,3],[1,-3]]
        b = [9,-3]

        c = [Fraction(x) for x in c]
        A = [[Fraction(x) for x in line] for line in A]
        b = [Fraction(x) for x in b]

        solution = [3,2]
        value = sum(solution[var] * c[var] for var in range(len(c)))

        self.assertEqual(solveStandardSimplex(c=c,A=A,b=b),(solution,value))

    def test_program3(self):
        c = [1,1]
        A = [[-1,2],[3,-1]]
        b = [4,5]

        c = [Fraction(x) for x in c]
        A = [[Fraction(x) for x in line] for line in A]
        b = [Fraction(x) for x in b]

        solution = [Fraction(14,5),Fraction(17,5)]
        value = sum(solution[var] * c[var] for var in range(len(c)))

        self.assertEqual(solveStandardSimplex(c=c,A=A,b=b),(solution,value))

    def test_program4(self):
        c = {'x':-1, 'y':1}
        A = [{'x':-1, 'y':3}, {'x':1, 'y':3}]
        b = [9,3]
        Aineq = ['<=', '>=']
        vartype = {'x':'n', 'y':'u'}

        solution = {'x':-3, 'y':2}
        value = 5

        self.assertEqual(solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=None, maximization=True, varclass=Fraction), (solution, value))

    def test_program5(self):
        names = ['x', 'y', 'z', 'w']
        c = dict(zip(names,[-1,-2,-3,4]))
        A = [dict(zip(names,[1,-1,3,-1])),dict(zip(names,[1,-2,1,-1])),dict(zip(names,[2,1,-7,-1]))]
        b = [5,2,3]

        Aineq = ['<=', '>=', '<=']
        vartype = dict(zip(names, ['p', 'p', 'p', 'n']))

        solution = dict(zip(names, [0, Fraction(13, 7), Fraction(4, 7), -Fraction(36, 7)]))
        value = sum(solution[var] * c[var] for var in names)

        self.assertEqual(solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=None, maximization=False, varclass=Fraction, pivotchoice=greedyStrategy), (solution, value))
        
    def test_program6(self):
        names = ['x', 'y', 'z']
        c = dict(zip(names,[1,1,3]))
        A = [dict(zip(names,[-1,2,0])),dict(zip(names,[3,-1,0])), dict(zip(names,[0,0,3]))]
        b = [4,5,9]

        Aineq = ['<=', '<=', '=']
        vartype = dict(zip(names, ['p', 'p', 'u']))

        solution = dict(zip(names,[Fraction(14,5),Fraction(17,5), 3]))
        value = sum(solution[var] * c[var] for var in names)

        self.assertEqual(solveNamedSimplex(c=c, A=A, b=b, Aineq=Aineq, vartype=vartype, c0=None, maximization=True, varclass=Fraction), (solution, value))

    def test_program7(self):
        c = [1,1]
        A = [[1,3],[1,-3]]
        b = [9,-3]

        c = [Fraction(x) for x in c]
        A = [[Fraction(x) for x in line] for line in A]
        b = [Fraction(x) for x in b]

        solution = [0,1]
        value = sum(solution[var] * c[var] for var in range(len(c)))

        self.assertEqual(solveStandardSimplex(c=c,A=A,b=b,cutpoint=Fraction(1,1)),(solution,value))

    def test_dot_product(self):
        a = [1,2,3,4]
        b = [-2,4,1,0]

        self.assertEqual(dot_product(a,b), 9)

    def test_injectPreProcess(self):
        A = [list(range(5))]
        base = [2,4]

        candidateBase = [1,0,4,1]
        order = [1,0,4]
        candidateBaseSet = set(order)

        self.assertEqual(injectPreProcess(base, candidateBase), (order, candidateBaseSet))

    def test_gramSchmidt(self):
        A = [[ 1, 1, 1, 4],
             [ 1,-1, 0, 1],
             [ 1, 1, 1, 0]]

        (_, A, _) = dictifySystem(None, A, None)

        order = [2, 1, 0, 3]
        base = [2, 1, 3]

        self.assertEqual(gramSchmidtFindRankIndices(A, order), base)


    def test_gramSchmidt2(self):
        A = [[ 1, 1, 1, 0],
             [ 1,-1, 0, 1],
             [ 1, 1, 1, 0]]

        (_, A, _) = dictifySystem(None, A, None)

        order = [3, 1, 0, 2]
        base = [3, 1]

        self.assertEqual(gramSchmidtFindRankIndices(A, order), base)

    def test_sequentialInject(self):
        c = [0, 0, 0, 1, -1]
        A = [[1, 0, 0, 2, -1/2],
             [0, 1, 0, 4, -2],
             [0, 0, 1, -8, 0]]
        b = [1, 4, 2]
        c0 = 1
        base = [0, 1, 2]
        cutpoint = None
        candidateBase = [4, 3, 2, 1]
        times = 1
        callback = dummy

        (c, A, b) = dictifySystem(c, A, b, removeZeros=True)

        resultc0 = sequentialInjectCandidateBase(c, A, b, c0, base, cutpoint, 
                                                 candidateBase, times, callback=callback)

        newA = [[0.5, 0.0, 0.0, 1.0, -0.25],
                [-2.0, 1.0, 0, 0, -1.0],
                [4.0, 0.0, 1.0, 0.0, -2.0]]
        newb = [0.5, 2.0, 6.0]
        newc = [-0.5, 0.0, 0.0, 0.0, -0.75]
        newc0 = 1.5
        newbase = [3, 1, 2]

        (newc, newA, newb) = dictifySystem(newc, newA, newb, removeZeros=True)

        self.assertEqual(newA, A)
        self.assertEqual(newb, b)
        self.assertEqual(newc, c)
        self.assertEqual(newc0, resultc0)
        self.assertEqual(newbase, base)

if __name__ == '__main__':
    unittest.main()
