from PMFpack import *
from pylab import *
from numpy import *
from scipy.special import erfinv, erf
from scipy.integrate import quad, odeint, ode
from time import time

class diffusion:
    """
    This class is designed as a test case for the various submodels in PMFpack.
    We solve an equation for the mixture fraction mean and second raw moment,
    given by the equations
    
        d<f>/dt = d^2<f>/dx^2
        
        d<f*f>/dt = d^2<f*f>/dx^2 - itau * ( <f*f> - <f>^2 )
    
    The domain is x = [-L/2, L/2] with Neumann boundary conditions at the ends. 
    The initial conditions are either 
        
        (init) <f> = 1-H(x), <f*f> = <f>
        
        or
               
        (init2) <f> = 1-H(x+L/4)+H(x-L/4), <f*f> = <f>
       
    where H(x) is the Heaviside function that is unity for positive x and zero otherwise. 
    init2 is periodic and can thus be solved with Fourier spectral methods. For stability
    the Heviside function is approxiamted with an error function (see init and init2).
        
    The central variance is given as sigma = <f*f> - <f>^2
    
    The parameter itau determines the rate of mixing compared to transport.
    Increase this parameter and note how the inhomogeneous modification to the CSD
    increases.     
    
    """
    def __init__(self,**kwargs):
        if(self.__class__.__name__ == 'diffusion'):
            self.par={
                        'N'   : 100,
                        'itau': 10.,
                        'L'   : 2.,
                        't'   : 0.1,
                        'Nt'  : 100
                            }
        self.par.update(kwargs)
        for name in self.par.keys():
            setattr(self, name, self.par[name])
        self.dy = self.L / self.N
        self.Np = self.N + 1
        self.Nm = self.N - 1
        self.Ntot = self.Np * 2
        self.itau = self.par['itau']
        self.tv = arange(0, self.t, self.t / self.Nt)
        self.xis = zeros(self.Ntot)
        self.xi = reshape(self.xis, (2, -1)) # xi[0]=fmean, xi[1]=second integer moment
        self.zz = linspace(-self.L / 2, self.L / 2, self.Np)
        self.wa = zeros(self.Np * 2) # work array
        
    def init(self,sigma=0.05):
        """
        Initialize odes
        """
        self.xi[0, :] = 1 - 0.5 * (1 + erf(self.zz / sigma))
        self.xi[1, :] = self.xi[0, :]
    
    def init2(self,sigma=0.05):
        """
        Initialize odes
        """
        NN = self.N / 2 + 1
        self.xi[0, :] = 1 - 0.5 * erf((self.zz + 0.5) / sigma) + 0.5 * erf((self.zz - 0.5) / sigma)
        self.xi[1, :] = self.xi[0,:]

    def func(self, xis, t):
        """
        Return right hand side of odes
        """
        self.wa = self.wa * 0.
        yv = reshape(self.wa, (2,-1))
        xi = reshape(xis, (2,-1))
        xi[:, 0] = 4. / 3. * xi[:, 1] - 1. / 3. * xi[:, 2]
        xi[:, -1] = 4. / 3. * xi[:, -2] - 1. / 3. * xi[:, -3]
        yv[0, 1:-1] = (xi[0, 2:] - 2. * xi[0, 1:-1] + xi[0, :-2]) / self.dy**2
        yv[1, 1:-1] = (xi[1, 2:] - 2. * xi[1, 1:-1] + xi[1, :-2]) / self.dy**2 - self.itau * (xi[1, 1:-1] - xi[0, 1:-1]**2)        
        return self.wa
    
    def solve(self):
        """
        Solve problem
        """
        self.z, self.infodict = odeint(self.func, self.xis, self.tv, full_output=True)
        sz = self.z.shape
        self.z = reshape(self.z, (sz[0], 2, self.Np))
        self.z[:, :, 0] = 4. / 3. * self.z[:, :, 1] - 1. / 3. * self.z[:, :, 2]
        self.z[:, :, -1] = 4. / 3. * self.z[:, :, -2] - 1. / 3. * self.z[:, :, -3]
        self.chi = self.itau * (self.z[:, 1, :] - self.z[:, 0, :]**2) # chi = 2N, i.e it includes the number 2 in the def
        
    def dxidy(self, t):
        """
        Return gradients of the solution, given as inlet at a certain time-step
        """
        dfdy = zeros(self.Np)
        dsdy = zeros(self.Np)
        dmdy = zeros(self.Np)
        dfdy[1:-1] = (self.z[t, 0, 2:] - self.z[t, 0, :-2]) / self.dy / 2.
        dmdy[1:-1] = (self.z[t, 1, 2:] - self.z[t, 1, :-2]) / self.dy / 2.
        dsdy[1:-1] = dmdy[1:-1] - 2. * self.z[t, 0, 1:-1] * dfdy[1:-1]
        return dfdy, dmdy, dsdy

def Conditionals(t, nx, sol, central=False, lookup=None, pmfclass=None):
    
    dfdy, dmdy, dsdy = sol.dxidy(t)
    fmean = sol.z[t, 0, nx]
    sigma = sol.z[t, 1, nx] - fmean**2
    if pmfclass == None:
        pmf = PMF(fmean, sigma, 0, 0, 0, 1)
    else:
        pmf = pmfclass
        pmf.set_parameters(fmean, sigma)
    pmf.chi = sol.chi[t, nx]
    pmf.set_fmean_gradient(dfdy[nx], 0, 0)
    pmf.set_sigma_gradient(dmdy[nx], 0, 0)

    if(lookup):
        lookuptable = GSLLookup(pmf.derivator)
        if (lookuptable.read_table(lookup) == -1):
            print "Generating lookup table ", lookup
            lookuptable.generate_table(50, 50, lookup)
        lookuptable(0, True) # Look up all values
        #lookuptable(0, True, True) # Look up all values using root-polishing
    else:
        pmf.compute(0, True)
    Nc = 100
    eta = linspace(0.01, 0.99, Nc)
    csd = zeros(Nc)
    pmf.CSD(eta, csd)
    hcsd = zeros(Nc)
    pmf.counterflow(eta, hcsd)
    pdf = zeros(Nc)
    pmf.PDF(eta, pdf)
    figure()
    plot(eta, pdf, 'k')
    plot(eta, hcsd, 'b')
    plot(eta, csd, 'r')
    legend(('PMF-PDF','PMF-CSD Homo','PMF-CSD Full'))
    return pmf
    
if __name__ == "__main__":
    sol = diffusion(t=0.5)
    sol.init()
    sol.solve()
    pmf = Conditionals(20, 20, sol, lookup="GSL_Lookuptable50.dat")
    #pmf = Conditionals(20, 20, sol)
    show()
    