from PMFpack import *
from numpy import linspace, zeros, sqrt, log, pi, exp
from pylab import plot, show, figure, legend, xlabel, ylabel, title, semilogy
from time import time
from scipy.stats import norm

pdf = norm().pdf
cdf = norm().cdf
Erfinv = norm().ppf
# Initialize a PMF class that will be used to compute conditional submodels
# The numbers here are taken from a simulation using mixinglayer.py
#fmean = 0.639756248702123
#sigma = 0.16674922174125 * fmean * (1 - fmean)

#fmean = 4.653900989e-03
fmean = 0.02
#sigma = 8.9500083540e-01 * fmean * (1 - fmean)
sigma = 0.999* fmean * (1 - fmean)

pmf = PMF(fmean, sigma, 1, 1, 0, 1)

# Get the necessary parameters from the same 1D simulation
pmf.set_fmean_gradient(-0.4178939524532, 0, 0)  # Last two are zero because the simulation is 1D
pmf.set_sigma_gradient(-0.5095137299633, 0, 0)
pmf.chi = 0.38430391476813  # Mean scalar dissipation rate (defined with the factor 2)
pmf.DT = 1.         # Turbulent diffusivity

# The simulations are run using integer moments
pmf.central = False

#pmf.reallocate_solver(pmf.froot, 1)
pmf.compute(0, True)

from scipy.special import erf, erfinv
def integrand(x, tau, alfa):
    X = 0.5*(1+erf((x - alfa) / (2*sqrt(tau))))
    return X * X * exp(- x * x / (2*(1-2*tau)));

def integrand2(x, tau, alfa):
    X = cdf((x*(sqrt(1-2*tau)) - alfa) / sqrt(2*tau))
    return X * X / sqrt(2*pi) * exp(- x * x / 2)

def integrand2X(x, tau, alfa):
    X = 0.5*(1+erf((x - alfa) / (2*sqrt(tau))))
    return X * X

def integrand3(x, tau, alfa):
    X = cdf((x*(sqrt(1-2*tau)) - alfa) / sqrt(2*tau))
    return (X-fmean) * (X-fmean) / sqrt(2*pi) * exp(- x * x / 2)

def integrand4(t, tau, alfa):
    x = t/(1-t*t)
    X = 0.5*(1+erf((x - alfa) / (2*sqrt(tau))))
    return X * X * exp(- x * x / (2*(1-2*tau)))*(1+t**2)/(1-t**2)**2;

def integrand5(t, tau, alfa):
    d = 1-t
    x = alfa + t / d
    X = 0.5*(1+erf(t / d / (2*sqrt(tau)))) / d
    return X * X * exp(- x * x / (2*(1-2*tau)));

def integrand6(t, tau, alfa):
    d = 1-t
    x = alfa - d / t
    X = 0.5*(1+erf(-d / t / (2*sqrt(tau)))) / t
    return X * X * exp(- x * x / (2*(1-2*tau)));

def dintegrand(t, tau, alfa):
    x = t/(1-t*t)
    tmp = (x - alfa) / (2*sqrt(tau))
    X = 0.5*(1+erf(tmp))
    return X * (alfa - x / (1-2*tau)) *exp(-tmp**2 - x * x / (2*(1-2*tau))) / sqrt(pi) / tau**(1.5) / sqrt(1-2*tau)*(1+t**2)/(1-t**2)**2

def lower(x, tau, alfa):
    return integrand(x, tau, alfa) - 1e-24

def lower2(x, tau, alfa):
    return integrand2(x, tau, alfa) - 1e-24
    
tau = pmf.tau
alfa = pmf.alfa
ff = sqrt(1-2*tau) 
ep = 2*sqrt(-(1-2*tau)*log(2*pi*(1-2*tau)*1e-32));
em = alfa+2*sqrt(tau)*erfinv(2e-16-1);
    
y = linspace(em, ep, 1000)
figure("original log ")
semilogy(y, integrand(y, tau, alfa))
semilogy(y, integrand2(y, tau, alfa), 'g')
semilogy(y, integrand2X(y, tau, alfa), 'r')

figure("original")
plot(y, integrand(y, tau, alfa))
plot(y, integrand2(y, tau, alfa), 'g')

#figure()
#plot(y, integrand2X(y, tau, alfa))

figure("transformed")
xx = linspace(-1, 1, 1000)
plot(xx[1:-1], integrand4(xx[1:-1], tau, alfa))

figure("transformed log")
semilogy(xx[1:-1], integrand4(xx[1:-1], tau, alfa))

figure("transformed 2")
x2 = linspace(0, 1, 1000)
plot(x2[:-1], integrand5(x2[:-1], tau, alfa))

figure("transformed 2 log")
semilogy(x2[:-1], integrand5(x2[:-1], tau, alfa))

figure("transformed 3")
x2 = linspace(0, 1, 1000)
plot(x2[1:], integrand6(x2[1:], tau, alfa))

figure("transformed 3 log")
semilogy(x2[1:], integrand6(x2[1:], tau, alfa))

y = linspace(em, ep, 1000)

figure("dfdtau original")
plot(xx[1:-1], dintegrand(xx[1:-1], tau, alfa))

show()


#from scipy.optimize import fsolve
#print fsolve(lower, -1, args=(tau, alfa))
#print fsolve(lower, 1, args=(tau, alfa))

#print fsolve(lower2, -1, args=(tau, alfa))
#print fsolve(lower2, 1, args=(tau, alfa))

#figure()
#y = linspace(-20, 20, 10000)
#semilogy(y, integrand(y, 0.499, alfa))
#semilogy(y, integrand2(y, 0.499, alfa), 'g')
##semilogy(y, integrand3(y, 0.0001, alfa), 'r')
##semilogy(y, integrand3(y, 0.4999, alfa), 'm')

#show()

##pmf.set_fdfsolver(1)  # Choose a Boost fdfsolver instead of GSL
#pmf.verbose = 0

## We can use a lookuptable for computing tau and its derivatives
#lookuptable = GSLLookup(pmf.derivator)
#t0 = time()
#Nt = 50
##lookuptable.generate_table(Nt, Nt, 'GSL_table_{0}.dat'.format(Nt)) # Create a small table for testing
##print "Generated Lookup table ", time() - t0
## When generated it can be read in here instead:
#lookuptable.read_table('GSL_table_50.dat')

## Test first by simply computing tau and its derivatives
#t0 = time()
#for i in range(100):      # 100 times for timing purpose
    #pmf.tau = 0.1         # Initial guess
    #pmf.compute(0, True)  # 0 for verbose, True for computing derivatives
#print "\nComputed tau  ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

### Find the same parameters using a lookuptable
#t0 = time()
#for i in range(100):
    #lookuptable(0, True)
#print "\nLookup without root-polishing:"
#print "Looked up tau ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

## Finally we test by looking up tau and then using a root-polishing routine for accuracy
## This basically means that the lookuptable is used to compute an initial guess
#t0 = time()
#for i in range(100):
    #lookuptable(0, True, True)
#print "\nLookup with root-polishing:"
#print "Looked up tau ", pmf.tau, " dtaudf ", pmf.dtaudf, " d2taudfdf ", pmf.d2taudfdf, " time ", (time() - t0) / 100

## When tau and its derivatives are computed we are able to
## compute a bunch of conditional models at practically no extra cost
#N = 1000
#eta = linspace(0.001, 0.999, N)
#csd = zeros(N)
#pmf.CSD(eta, csd)
#cf = zeros(N)
#pmf.counterflow(eta, cf)
#pdf = zeros(N)
#pmf.PDF(eta, pdf)

#plot(eta, csd, 'b', eta, cf, 'r', eta, pdf, 'k')
#title("CSD and PDF")
#xlabel("$\eta$", fontsize='large')
#legend(('PMF-CSD Homogen (Counterflow)','PMF-CSD Inhomogeneous', 'PMF-PDF'))
#cv = zeros((3, N))
#pmf.CV(eta, cv)
#figure()
#plot(eta, cv[0, :])
#xlabel("$\eta$", fontsize='large')
#title("Conditional mean velocity")
##show()

#tau = zeros((lookuptable.Nf, lookuptable.Ns))
#dtaudf = zeros((lookuptable.Nf, lookuptable.Ns))
#dtauds = zeros((lookuptable.Nf, lookuptable.Ns))
#dtaudfdf = zeros((lookuptable.Nf, lookuptable.Ns))
#dtaudsds = zeros((lookuptable.Nf, lookuptable.Ns))
#tau_table = lookuptable.get_table(0)
#dtaudf_table = lookuptable.get_table(1)
#dtauds_table = lookuptable.get_table(2)
#dtaudfdf_table = lookuptable.get_table(3)
#dtaudsds_table = lookuptable.get_table(4)
#for i in range(lookuptable.Nf):
    #for j in range(lookuptable.Ns):
        #tau[i, j] = double_getitem(doublep_getitem(tau_table, i), j)
        #dtaudf[i, j] = double_getitem(doublep_getitem(dtaudf_table, i), j)
        #dtauds[i, j] = double_getitem(doublep_getitem(dtauds_table, i), j)
        #dtaudfdf[i, j] = double_getitem(doublep_getitem(dtaudfdf_table, i), j)
        #dtaudsds[i, j] = double_getitem(doublep_getitem(dtaudsds_table, i), j)

#fm = zeros(lookuptable.Nf)
#for i in range(lookuptable.Nf):
    #fm[i] = double_getitem(lookuptable.fm, i)
    
#Is = zeros(lookuptable.Ns)
#for i in range(lookuptable.Ns):
    #Is[i] = double_getitem(lookuptable.Is, i)

#from pylab import contourf, colorbar
#figure('tau')
#contourf(fm, Is, tau, 100)
#colorbar()

#figure('dtaudf')
#contourf(fm, Is, dtaudf, 100)
#colorbar()

#figure('dtau_ds')
#contourf(fm, Is, dtauds, 100)
#colorbar()

#figure('dtau_dfdf')
#contourf(fm, Is, dtaudfdf, 100)
#colorbar()

#figure('dtau_dsds')
#contourf(fm, Is, dtaudsds, 100)
#colorbar()

##show()
