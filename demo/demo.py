from PMFpack import *
from numpy import linspace, zeros
from pylab import plot, show
from time import time

fm = 0.70837758
sigma = 0.194360748 * fm * (1 - fm)
pmf = PMF(fm, sigma, 0, 0, 0, 1)
pmf.set_fmean_gradient(-0.3496497, 0, 0)
pmf.set_sigma_gradient(-0.4294126, 0, 0)
pmf.chi = 1.003770
pmf.DT = 1.
pmf.central = False
#pmf.reallocate_solver(pmf.froot, 1)
pmf.compute(0, True)
pmf.tau = 0.1

#pmf.set_fdfsolver(1)
pmf.verbose = 0

# Test lookuptable
lookuptable = GSLLookup(pmf.derivator)
t0 = time()
lookuptable.generate_table(10, 10, 'GSL_table.dat') # Create a small table for testing
print "Generated Lookup table ", time() - t0

t0 = time()
for i in range(100): 
    pmf.tau = 0.1
    pmf.compute(0, True)
print "\nComputed tau  ", pmf.tau, " time ", (time() - t0) / 100

# When generated it can be read in here instead:
#lookuptable.read_table('GSL_table.dat')
t0 = time()
for i in range(100):
    lookuptable(0, True)
print "\nLookup without root-polishing:"
print "Looked up tau ", pmf.tau, " time ", (time() - t0) / 100

t0 = time()
for i in range(100):
    lookuptable(0, True, True)
print "\nLookup with root-polishing:"
print "Looked up tau ", pmf.tau, " time ", (time() - t0) / 100
# Test some conditional models
N = 1000
eta = linspace(0.001, 0.999, N)
csd = zeros(N)
pmf.CSD(eta, csd)
cf = zeros(N)
pmf.counterflow(eta, cf)

#plot(eta, csd, 'b', eta, cf, 'r')
#show()
cv = zeros((3, N))
pmf.CV(eta, cv)

