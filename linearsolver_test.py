import unittest

from linearsolver import *
import simplex

class TestLinearsolver(unittest.TestCase):

    def test_primal(self):
        (x,v) = solveBirkhoffPrimal(2, 15, fractions.Fraction(17,10), callback=simplex.generateModPrinter(10))
        self.assertEqual(v, fractions.Fraction(-693546131200007083960247889737617216097620367281230353754881253059537327, 899499953539818762045703020105883975680000000000000000000000000000000000))

    def test_dual(self):
        (x,v) = solveBirkhoffDual(2, 15, fractions.Fraction(17,10), callback=simplex.generateModPrinter(10))

        self.assertEqual(v, fractions.Fraction(-693546131200007083960247889737617216097620367281230353754881253059537327, 899499953539818762045703020105883975680000000000000000000000000000000000))

    def test_primal_vs_dual(self):
        (xp, vp) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10), callback=simplex.generateModPrinter(10))
        (xd, vd) = solveBirkhoffDual(2, 17, fractions.Fraction(17,10), callback=simplex.generateModPrinter(10))
        self.assertEqual(vp, vd)

    def test_primal_vs_dual2(self):
        (xp, vp) = solveBirkhoffPrimal(4, 15, fractions.Fraction(16,10),
                                       callback=simplex.generateModPrinter(10),
                                       lset=set([0,4]))
        (xd, vd) = solveBirkhoffDual(4, 15, fractions.Fraction(16,10),
                                     callback=simplex.generateModPrinter(10),
                                     lset=set([0,4]))
        self.assertEqual(vp, vd)
        self.assertTrue(vp > 0)

    def test_primal_vs_dual3(self):
        (x, v) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                     callback=simplex.generateModPrinter(10),
                                     kostkaRounder=dummy)
        (xp, vp) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                       callback=simplex.generateModPrinter(10),
                                       kostkaRounder=generateMildKostkaRounder(fractions.Fraction(1,100)))
        (xd, vd) = solveBirkhoffDual(2, 17, fractions.Fraction(17,10),
                                     callback=simplex.generateModPrinter(10),
                                     kostkaRounder=generateMildKostkaRounder(fractions.Fraction(1,100)))
        (xp2, vp2) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                         callback=simplex.generateModPrinter(10),
                                         kostkaRounder=generateAggressiveKostkaRounder(20))
        (xd2, vd2) = solveBirkhoffDual(2, 17, fractions.Fraction(17,10),
                                       callback=simplex.generateModPrinter(10),
                                       kostkaRounder=generateAggressiveKostkaRounder(20))
        self.assertEqual(vp, vd)
        self.assertEqual(vp2, vd2)
        self.assertTrue(vp <= v)
        self.assertTrue(vp2 <= v)

    def test_primal_vs_dual4(self):
        (x, v) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                     callback=simplex.generateModPrinter(10),
                                     negativeBelly=False,
                                     negativeLeg=False)
        (xp, vp) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                       callback=simplex.generateModPrinter(10),
                                       negativeBelly=True,
                                       negativeLeg=False)
        (xd, vd) = solveBirkhoffDual(2, 17, fractions.Fraction(17,10),
                                     callback=simplex.generateModPrinter(10),
                                     negativeBelly=True,
                                     negativeLeg=False)
        (xp2, vp2) = solveBirkhoffPrimal(2, 17, fractions.Fraction(17,10),
                                         callback=simplex.generateModPrinter(10),
                                         negativeBelly=False,
                                         negativeLeg=True)
        (xd2, vd2) = solveBirkhoffDual(2, 17, fractions.Fraction(17,10),
                                       callback=simplex.generateModPrinter(10),
                                       negativeBelly=False,
                                       negativeLeg=True)
        self.assertEqual(vp, vd)
        self.assertEqual(vp2, vd2)
        self.assertTrue(vp <= v)
        self.assertTrue(vp2 <= v)

    def test_fragdual(self):
        (problem, sold, defragsold) = solveFragBirkhoffDual(2, 15, fractions.Fraction(17,10),
                                                       callback=simplex.generateModPrinter(10), save=False)

        feasible = defragsold['feasible']
        (x,v) = defragsold['solution']

        self.assertEqual(v, fractions.Fraction(-693546131200007083960247889737617216097620367281230353754881253059537327, 899499953539818762045703020105883975680000000000000000000000000000000000))

        self.assertEqual(feasible, True)

    def test_primal_vs_fragdual(self):
        (xp, vp) = solveBirkhoffPrimal(4, 15, fractions.Fraction(16,10),
                                       callback=simplex.generateModPrinter(10),
                                       lset=set([0,4]))
        (problem, sold, defragsold) = solveFragBirkhoffDual(4, 15, fractions.Fraction(16,10),
                                                         callback=simplex.generateModPrinter(10),
                                                         lset=set([0,4]), save=False)
        feasible = defragsold['feasible']
        (xd,vd) = defragsold['solution']
        self.assertEqual(vp, vd)
        self.assertTrue(vp > 0)
        self.assertEqual(feasible, True)

if __name__ == '__main__':
    unittest.main()
