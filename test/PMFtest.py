from PMFpack import *

pmf = PMF(0.5, 0.1, 0, 0)
pmf.initial_guess(0.1)
pmf.verbose = 1;
pmf.reallocate_solver(pmf.froot, 1)
pmf.computeall()