##############
#  Standard  #
##############

parameters = [( 0, 19, fractions.Fraction(149,100))]
parameters = [( 2, 29, fractions.Fraction(169,100))]
parameters = [( 4, 29, fractions.Fraction(172,100))]
parameters = [( 6, 39, fractions.Fraction(178,100))]
parameters = [( 8, 39, fractions.Fraction(180,100))]
parameters = [(10, 49, fractions.Fraction(182,100))]
parameters = [(12, 59, fractions.Fraction(185,100))]
parameters = [(14, 79, fractions.Fraction(187,100))]

def runPrimal():
    for (lmax, k0, c) in parameters:
        print('\n\nStandard\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveBirkhoffPrimalAndSave(lmax=lmax, k0=k0, c=c, cutpoint=0, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=range(2,2*(lmax+k0),2), fileprefix='even_')

def runDual():
    for (lmax, k0, c) in parameters:
        print('\n\nStandard\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveBirkhoffDualAndSave(lmax=lmax, k0=k0, c=c, cutpoint=None, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=range(2,2*(lmax+k0),2), fileprefix='even_')


##########################
#  Frag/defrag Standard  #
##########################

parameters = [( 0,  19, fractions.Fraction(149,100))]
parameters = [( 2,  29, fractions.Fraction(169,100))]
parameters = [( 4,  29, fractions.Fraction(172,100))]
parameters = [( 6,  39, fractions.Fraction(178,100))]
parameters = [( 8,  39, fractions.Fraction(180,100))]
parameters = [(10,  49, fractions.Fraction(182,100))]
parameters = [(12,  59, fractions.Fraction(185,100))]
parameters = [(14,  79, fractions.Fraction(187,100))]
parameters = [(16,  99, fractions.Fraction(189,100))]
parameters = [(18, 109, fractions.Fraction(190,100))]
parameters = [(20,  89, fractions.Fraction(190,100))]
parameters = [(22, 109, fractions.Fraction(191,100))]
parameters = [(24, 129, fractions.Fraction(192,100))]
parameters = [(26, 119, fractions.Fraction(192,100))]
parameters = [(28, 149, fractions.Fraction(193,100))]
parameters = [(30, 139, fractions.Fraction(193,100))]
parameters = [(32, 189, fractions.Fraction(194,100))]

def runFragDual():
    for (lmax, k0, c) in parameters:
        print('\n\nFrag/defrag Standard\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveFragBirkhoffDual(lmax=lmax, k0=k0, c=c, cutpoint=None, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=range(2,2*(lmax+k0),2), save=True, lightweight=1, verbose=True, verbosePeriod=100, fileprefix='even_')


#########################################################################
#  Frag/defrag Joint Large Leg Mild 0.01 noRescale Clever Skip 4-30-50  #
#########################################################################
    
parameters = [( 2,  29, fractions.Fraction(148,100))]
parameters = [( 4,  29, fractions.Fraction(159,100))]
parameters = [( 6,  29, fractions.Fraction(171,100))]
parameters = [( 8,  39, fractions.Fraction(176,100))]
parameters = [(10,  49, fractions.Fraction(180,100))]
parameters = [(12,  49, fractions.Fraction(183,100))]
parameters = [(14,  79, fractions.Fraction(186,100))]
parameters = [(16,  99, fractions.Fraction(188,100))]
parameters = [(18,  89, fractions.Fraction(189,100))]
parameters = [(20,  99, fractions.Fraction(190,100))]
parameters = [(22, 119, fractions.Fraction(191,100))]
parameters = [(24, 109, fractions.Fraction(191,100))]
parameters = [(26, 129, fractions.Fraction(192,100))]
parameters = [(28, 179, fractions.Fraction(193,100))]
parameters = [(30, 149, fractions.Fraction(193,100))]
parameters = [(32, 139, fractions.Fraction(193,100))]
parameters = [(34, 189, fractions.Fraction(194,100))]
parameters = [(36, 169, fractions.Fraction(194,100))]
parameters = [(38, 169, fractions.Fraction(194,100))]
parameters = [(40, 169, fractions.Fraction(194,100))]
parameters = [(42, 229, fractions.Fraction(195,100))]
parameters = [(44, 219, fractions.Fraction(195,100))]
parameters = [(46, 219, fractions.Fraction(195,100))]
parameters = [(48, 209, fractions.Fraction(195,100))]
parameters = [(50, 199, fractions.Fraction(195,100))]
parameters = [(52, 199, fractions.Fraction(195,100))]
parameters = [(54, 309, fractions.Fraction(196,100))]
parameters = [(56, 289, fractions.Fraction(196,100))]
parameters = [(58, 279, fractions.Fraction(196,100))]
parameters = [(60, 269, fractions.Fraction(196,100))]
parameters = [(62, 269, fractions.Fraction(196,100))]
parameters = [(64, 269, fractions.Fraction(196,100))]
parameters = [(66, 259, fractions.Fraction(196,100))]
parameters = [(68, 259, fractions.Fraction(196,100))]
parameters = [(70, 539, fractions.Fraction(197,100))]
parameters = [(70, 539, fractions.Fraction(19701,10000))]
parameters = [(74, 469, fractions.Fraction(1971,1000))]

def runFragDual():
    cskipify(parameters, firstinc=4, smallskip=30, bigskip=50)
    for (lmax, k0, c, mset) in parameters:
        print('\n\nFrag/defrag Joint Large Leg Mild 0.01 no rescale Clever Skip 4-30-50\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveFragBirkhoffDual(lmax=lmax, k0=k0, c=c, cutpoint=None, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=mset, kostkaRounder=generateMildKostkaRounder(fractions.Fraction(1,100), rescale=False), jointLargeLeg=True, save=True, lightweight=1, verbose=True, verbosePeriod=100, fileprefix='even_jointLargeLeg_mild0.01noRescale_cleverSkip4-30-50_', num_workers=12)

def runFragDualLightweight():
    cskipify(parameters, firstinc=4, smallskip=30, bigskip=50)
    for (lmax, k0, c, mset) in parameters:
        print('\n\nFrag/defrag Joint Large Leg Mild 0.01 no rescale Clever Skip 4-30-50\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveFragBirkhoffDual(lmax=lmax, k0=k0, c=c, cutpoint=None, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=mset, kostkaRounder=generateMildKostkaRounder(fractions.Fraction(1,100), rescale=False), jointLargeLeg=True, save=True, lightweight=2, verbose=True, verbosePeriod=100, fileprefix='even_jointLargeLeg_mild0.01noRescale_cleverSkip4-30-50_', num_workers=12)

def runFragFirstPhase():
    cskipify(parameters, firstinc=4, smallskip=30, bigskip=50)
    for (lmax, k0, c, mset) in parameters:
        print('\n\nFrag/defrag Joint Large Leg Mild 0.01 no rescale Clever Skip 4-30-50 (First phase only)\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
        sol = solveFragBirkhoffDual(lmax=lmax, k0=k0, c=c, cutpoint=None, callback=simplex.joinFunctions(simplex.generateModPrinter(),simplex.generateAverageSizeTracker(),simplex.generateTotalSizeTracker()), pivotchoice=simplex.greedyStrategy, mset=mset, kostkaRounder=generateMildKostkaRounder(fractions.Fraction(1,100), rescale=False), jointLargeLeg=True, save=True, lightweight=2, verbose=True, verbosePeriod=100, fileprefix='even_jointLargeLeg_mild0.01noRescale_cleverSkip4-30-50_', num_workers=12, runFragmentedPhaseOnly=True)

def runFragSecondPhase(filename):
    sold = load(filename)
    resumeFragBirkhoffDual(sold)
