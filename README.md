
# General 

This repository contains the linear programming package used to obtain
the results of the paper "Tighter Bounds on the Independence Number of
the Birkhoff Graph" by Leonardo Nagami Coregliano and Fernando Granha
Jeronimo. Please, see LICENSE file for legal information.

# File Structure 

- CompactSet.py: meta-data compression module.
- LICENSE: legal information 
- linearsolver.py: interface module containing linear program generators and solver wrappers.
- linearsolver_test.py: linear solver test module.
- Partition.py: Partition class and associated functionalities.
- partition_test.py: partition test module.
- run.py: file containing example of how to invoke the linear solver (see more details in the next section).
- simplex.py: simplex linear programming solver with support to rational computation.
- simplex_test.py: simplex solver test module.

# Running Instructions

Here is an example on how to use the code snipets from run.py to
execute the standard primal program with parameters l_0=4, k=29,
c=1.72.

```python
>>> from linearsolver import *
>>> parameters = [( 4, 29, fractions.Fraction(172,100))]
>>> def runPrimal():
...    for (lmax, k0, c) in parameters:
...        print('\n\nStandard\nlmax = %d, k0 = %d, c = %.21g' % (lmax, k0, c))
...        sol = solveBirkhoffPrimalAndSave(lmax=lmax, k0=k0, c=c, cutpoint=0, callback=simplex.generateModPrinter(), pivotchoice=simplex.greedyStrategy, mset=range(2,2*(lmax+k0),2), fileprefix='even_')
>> runPrimal()
```

The run.py file has three sections, each with its own functions:
  - Stardard: original primal/dual without speed up heuristics.
    - runPrimal(): creates, solves and saves the primal problem.
    - runDual(): creates, solves and saves the dual problem.
  - Frag/defrag Standard: dual with only fragment heuristic (see paper).
    - runFragDual(): creates, solves and saves the dual fragmented problem.
  - Frag/defrag Joint Large Leg Mild 0.01 noRescale Clever Skip 4-30-50:
     dual with all heuristics (see paper).
    - runFragDual(): creates, solves and saves the dual problem with all speed
                     up heuristics, producing meta-data for later analysis
		     (not recommended for large values of lmax=l_0).
    - runFragDualLightweight(): creates, solves and saves the dual problem with all speed
                                up heuristics, producing less meta-data for later analysis.
    - runFragFirstPhase(): creates, solves and saves the first phase (the fragment problem) of the dual
                           problem. It generates a frag solution file to be used in the second phase.   
    - runFragSecondPhase(filename): checks if dual frag solution file from the first phase satifies
                                    all the original linear programming restrictions.

Warning: the file run.py is not intended to be executed directly but to have snipets of its code
         to be copied and executed.

More specific details about the parameters of the module functions can
be found using the usual python's help fuction.

# Solution files

The data folder contains solution files outputted by the program. The solution files are subclassified into the folders Standard, FragHeuristic and AllHeuristics corresponding to the three types of examples in run.py.

To inspect a particular solution file, you can use the functions load, loadall, printsol, printrawsol, printallsol, printallrawsol from linearsolver.py (for more details on these functions, use python's help function). The variables of the problem are described in the article.

Example:
```python
>>> from linearsolver import *
>>> sol = load('data/Standard/even_dualsol_14_79_1.87_mset.pkl')
>>> printsol(sol)
```
