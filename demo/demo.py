from PMFpack import *
from numpy import linspace, zeros
from pylab import plot, show, figure, legend, xlabel, ylabel, title
from time import time

# Initialize a PMF class that will be used to compute conditional submodels
# The numbers here are taken from a simulation using mixinglayer.py
fmean = 0.70837758
sigma = 0.194360748 * fmean * (1 - fmean)
pmf = PMF(fmean, sigma, 0, 0, 0, 2)

# Get the necessary parameters from the same 1D simulation
pmf.set_fmean_gradient(-0.3496497, 0, 0)  # Last two are zero because the simulation is 1D
pmf.set_sigma_gradient(-0.4294126, 0, 0)
pmf.chi = 1.003770  # Mean scalar dissipation rate (defined with the factor 2)
pmf.DT = 1.         # Turbulent diffusivity

# The simulations are run using integer moments
pmf.central = False

#pmf.reallocate_solver(pmf.froot, 1)
pmf.compute(0, True)

#pmf.set_fdfsolver(1)  # Choose a Boost fdfsolver instead of GSL
pmf.verbose = 0

# We can use a lookuptable for computing tau and its derivatives
lookuptable = GSLLookup(pmf.derivator)
t0 = time()
Nt = 10
lookuptable.generate_table(Nt, Nt, 'GSL_table_{0}.dat'.format(Nt)) # Create a small table for testing
print "Generated Lookup table ", time() - t0
# When generated it can be read in here instead:
#lookuptable.read_table('GSL_table_10.dat')

# Test first by simply computing tau and its derivatives
t0 = time()
for i in range(1):      # 100 times for timing purpose
    pmf.tau = 0.1         # Initial guess
    pmf.compute(0, True)  # 0 for verbose, True for computing derivatives
print "\nComputed tau  ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

## Find the same parameters using a lookuptable
t0 = time()
for i in range(1):
    lookuptable(0, True)
print "\nLookup without root-polishing:"
print "Looked up tau ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

# Finally we test by looking up tau and then using a root-polishing routine for accuracy
# This basically means that the lookuptable is used to compute an initial guess
t0 = time()
for i in range(1):
    lookuptable(0, True, True)
print "\nLookup with root-polishing:"
print "Looked up tau ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

# When tau and its derivatives are computed we are able to
# compute a bunch of conditional models at practically no extra cost
N = 1000
eta = linspace(0.001, 0.999, N)
csd = zeros(N)
pmf.CSD(eta, csd)
cf = zeros(N)
pmf.counterflow(eta, cf)
pdf = zeros(N)
pmf.PDF(eta, pdf)

plot(eta, csd, 'b', eta, cf, 'r', eta, pdf, 'k')
title("CSD and PDF")
xlabel("$\eta$", fontsize='large')
legend(('PMF-CSD Homogen (Counterflow)','PMF-CSD Inhomogeneous', 'PMF-PDF'))
cv = zeros((3, N))
pmf.CV(eta, cv)
figure()
plot(eta, cv[0, :])
xlabel("$\eta$", fontsize='large')
title("Conditional mean velosity")
show()
