from PMFpack import *
from numpy import linspace, zeros
from pylab import plot, show
fm = 0.70837758
sigma = 0.194360748 * fm * (1 - fm)
pmf = PMF(fm, sigma, 0, 0, 0, 1)
pmf.set_fmean_gradient(-0.3496497, 0, 0)
pmf.set_sigma_gradient(-0.4294126, 0, 0)
pmf.chi = 1.003770
pmf.DT = 1.
pmf.central = False

#double_setitem(pmf.tau, 3, 0.1) # Initial guess for tau
pmf.verbose = 1
#pmf.reallocate_solver(pmf.froot, 1)
pmf.compute(0, True)

N = 1000
eta = linspace(0.001, 0.999, N)
csd = zeros(N)
pmf.CSD(eta, csd)
cf = zeros(N)
pmf.counterflow(eta, cf)

plot(eta, csd, 'b', eta, cf, 'r')
show()
cv = zeros((3, N))
pmf.CV(eta, cv)