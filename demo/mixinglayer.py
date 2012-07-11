from PMFpack import *
from pylab import *
from numpy import *
from scipy.special import erfinv, erf
from scipy.integrate import quad, odeint, ode
from scipy.optimize import fsolve
from time import time

class mixinglayer:
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
    def __init__(self, **kwargs):
        if self.__class__.__name__ == 'mixinglayer':
            self.par={
                        'N'   : 100,
                        'itau': 10.,
                        'L'   : 8.,
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
        self.xi[0, :] = (1 - 0.5 * (1 + erf(self.zz / sigma))) * 0.98 + 0.01
        self.xi[1, :] = self.xi[0, :]
    
    def init2(self,sigma=0.05):
        """
        Initialize odes
        """
        NN = self.N / 2 + 1
        self.xi[0, :] = 1 - 0.5 * erf((self.zz + 0.5) / sigma) + 0.5 * erf((self.zz - 0.5) / sigma)
        self.xi[1, :] = self.xi[0, :]**1.01

    def func(self, xis, t):
        """
        Return right hand side of odes
        """
        wa = self.wa[:2*self.Np] * 0.
        yv = reshape(wa, (2, -1))
        xi = reshape(xis, (2, -1))
        xi[:, 0] = xi[:, 1]
        xi[:, -1] = xi[:, -2]
        #xi[:, 0] = 4. / 3. * xi[:, 1] - 1. / 3. * xi[:, 2]
        #xi[:, -1] = 4. / 3. * xi[:, -2] - 1. / 3. * xi[:, -3]
        #yv[0, 0] = (xi[0, 1] - xi[0, 0]) / self.dy
        #yv[1, 0] = (xi[1, 1] - xi[1, 0]) / self.dy - self.itau * (xi[1, 0] - xi[0, 0]**2)        
        #yv[0, -1] = -(xi[0, -1] - xi[0, -2]) / self.dy
        #yv[1, -1] = -(xi[1, -1] - xi[1, -2]) / self.dy - self.itau * (xi[1, -1] - xi[0, -1]**2)
        yv[0, 1:-1] = (xi[0, 2:] - 2. * xi[0, 1:-1] + xi[0, :-2]) / self.dy**2
        yv[1, 1:-1] = (xi[1, 2:] - 2. * xi[1, 1:-1] + xi[1, :-2]) / self.dy**2 - self.itau * (xi[1, 1:-1] - xi[0, 1:-1]**2)
        return wa
    
    def solve(self):
        """
        Solve problem
        """
        self.z, self.infodict = odeint(self.func, self.xis, self.tv, full_output=True)
        sz = self.z.shape
        self.z = reshape(self.z, (sz[0], 2, self.Np))
        #self.z[:, :, 0] = self.z[:, :, 1] 
        #self.z[:, :, -1] = self.z[:, :, -2]
        #self.z[:, :, 0] = 4. / 3. * self.z[:, :, 1] - 1. / 3. * self.z[:, :, 2]
        #self.z[:, :, -1] = 4. / 3. * self.z[:, :, -2] - 1. / 3. * self.z[:, :, -3]
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
        

def Conditionals(t, nx, sol, central=False, lookup=None, pmfclass=None, _plot=False):
    """Test of conditional models."""
    # Get the solution from mixinglayer at a certain time t and position nx
    dfdy, dmdy, dsdy = sol.dxidy(t)
    fmean = sol.z[t, 0, nx]
    sigma = sol.z[t, 1, nx] - fmean**2
    
    if fmean < 1e-6 or fmean > 0.999999:
        return None
    xx = sigma / fmean / (1 - fmean)
    if xx < 1e-8 or xx > 0.999999:
        return None
    
    # Create a PMF class and set all parameters according to simulation
    if pmfclass == None:
        pmf = PMF(fmean, sigma, 0, 0, 0, 2)
    else:
        pmf = pmfclass
        pmf.set_parameters(fmean, sigma)
    pmf.chi = sol.chi[t, nx]
    pmf.set_fmean_gradient(dfdy[nx], 0, 0)
    pmf.set_sigma_gradient(dmdy[nx], 0, 0)

    # Solve for tau and its derivatives
    if(lookup):
        lookuptable = GSLLookup(pmf.derivator)
        if (lookuptable.read_table(lookup) == -1):
            print "Generating lookup table ", lookup
            lookuptable.generate_table(50, 50, lookup)
        lookuptable(0, _plot) # Look up all values
        #lookuptable(0, True, True) # Look up all values using root-polishing
    else:
        pmf.compute(0, True)
        
    if _plot:
        # Compute and plot conditional models
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
    
store = {}
class reactive_mixinglayer(mixinglayer):
    def __init__(self, **kwargs):
        if self.__class__.__name__ == 'reactive_mixinglayer':
            self.par={
                        'N'   : 100,
                        'vode': True,
                        'itau': 10.,
                        'L'   : 2.,
                        't'   : 0.1,
                        'Nt'  : 100,
                        'r'   : 3.,
                        'alpha': 0.87,
                        'Ze'   : 4.0,
                        'A'    : 5e7,
                        'kappa': 0.01,
                        'N_CMC': 20}
        mixinglayer.__init__(self, **kwargs)
        self.Ntot = self.Np * (2 + self.N_CMC)
        self.eta = linspace(0, 1, self.N_CMC)
        self.deta = self.eta[1] - self.eta[0]
        self.csd = zeros(self.N_CMC - 2)
        self.cv  = zeros((3, self.N_CMC-2))
        self.pmf = [PMF(0.5, 0.1, 0, 0, 0, 0) for i in range(self.Np)]
        self.xis = zeros(self.Ntot)
        self.xi = reshape(self.xis, (self.N_CMC + 2, -1)) # xi[0] = fmean, xi[1] = im, xi[2] = T(eta=0), xi[3] = T(eta=deta)
        self.wa = zeros(self.Np * (self.N_CMC + 2)) # work array
        
    def init(self, sigma=0.05):
        mixinglayer.init(self, sigma)
        # init CMC solution to equilibrium
        p = self.par
        if p['r'] == 1:
            a = p['r'] / (p['r'] + 1.)**2 - p['kappa']
            b = (self.eta / (p['r'] + 1.) - 1. / (p['r'] + 1.) - self.eta * p['r'] / (p['r'] + 1.))
            c = self.eta - self.eta**2
            T = (-b - sqrt(b**2 - 4 * a * c)) / 2 / a
        else:
            T = fsolve(self.q0, self.eta)
        for i in range(self.Np):
            self.xi[2:, i] = T[:]

    def q0(self, x):
        r, kappa = self.par['r'], self.par['kappa']
        return (self.eta - x / (r+1)) * (1 - self.eta - r * x / (r + 1))**r - kappa * x**(r+1)

    def funcT(self, xis, t):
        """
        Return right hand side of odes
        """
        self.wa = self.wa * 0.
        yv = reshape(self.wa, (self.N_CMC + 2, -1))
        xi = reshape(xis, (self.N_CMC + 2, -1))
        # All equation Neumann at the ends
        #xi[:, 0] = 4. / 3. * xi[:, 1] - 1. / 3. * xi[:, 2]
        #xi[:, -1] = 4. / 3. * xi[:, -2] - 1. / 3. * xi[:, -3]
        xi[:, 0] = xi[:, 1]
        xi[:, -1] = xi[:, -2]
        # CMC has homogeneous Dirichlet at ends
        xi[2, :] = 0
        xi[-1, :] = 0
        # Mean and variance
        yv[0, 1:-1] = (xi[0, 2:] - 2. * xi[0, 1:-1] + xi[0, :-2]) / self.dy**2
        yv[1, 1:-1] = (xi[1, 2:] - 2. * xi[1, 1:-1] + xi[1, :-2]) / self.dy**2 - self.itau * (1. - exp(-5 * t)) * (xi[1, 1:-1] - xi[0, 1:-1]**2)
        # CMC
        yv[3:(self.N_CMC+1), 1:-1] = (xi[3:(self.N_CMC+1), 2:] - 2. * xi[3:(self.N_CMC+1), 1:-1] + xi[3:(self.N_CMC+1), :-2]) / self.dy**2

        for j in range(1, self.Np-1):            
            dfdy = (xi[0, j+1] - xi[0, j-1]) / self.dy / 2.
            dmdy = (xi[1, j+1] - xi[1, j-1]) / self.dy / 2.
            f, im = xi[0:2, j]
            s = im - f * f
            Is = s / f / (1 - f)
            if f > 0.01 and f < 0.99 and Is > 0.01 and Is < 0.99: 
                self.pmf[j].set_parameters(f, s)
                self.pmf[j].set_fmean_gradient(dfdy, 0, 0)
                self.pmf[j].set_sigma_gradient(dmdy, 0, 0)
                self.pmf[j].chi = self.itau * (1. - exp(-5 * t)) * s
                self.pmf[j].compute(0, True)
                #print 'Computing PMF'
        Ze = self.par['Ze']
        r = self.par['r']
        A = self.par['A']
        kappa = self.par['kappa']
        alfa = self.par['alpha']
        for j in range(1, self.Np-1):            
            f, im = xi[0:2, j]
            s = im - f * f
            Is = s / f / (1 - f)
            if f > 0.01 and f < 0.99 and Is > 0.01 and Is < 0.99: 
                self.csd[:] = self.pmf[j].chi
                self.pmf[j].counterflow(self.eta[1:-1], self.csd)
                self.pmf[j].CV(self.eta[1:-1], self.cv)
                rx = 1.
            else:
                self.csd[:] = 0.
                self.cv[0, :] = 0.
                rx = 0.
            yv[3:-1, j] += self.csd[:] / 2. * (xi[4:, j] - 2 * xi[3:-1, j] + xi[2:-2, j]) / self.deta**2
            yv[3:-1, j] += rx * ( (r + 1) * A * exp(- Ze / alfa - Ze * (1 - xi[3:-1, j]) / (1 - alfa * (1 - xi[3:-1, j])))
                * ((self.eta[1:-1] - xi[3:-1, j] / (r + 1)) * (1 - self.eta[1:-1] - r * xi[3:-1, j] / (r+1))**r - kappa * xi[3:-1, j]**(r+1)) )
            yv[3:-1, j] -= 2. * self.cv[0, :] * (xi[3:-1, j+1] - xi[3:-1, j-1]) / 2. / self.dy
        print t
        return self.wa
        
    # ode requires the parameters in opposite order from odeint
    def funcT2(self, t, xis):
        return self.funcT(xis, t)
        
    def solve(self):
        """
        Solve problem
        """
        if self.par['vode']:
            ode_runner = ode(self.funcT2)
            ode_runner.set_integrator('vode', method='bdf', rtol=1e-6, atol=1e-6, with_jacobian=False, max_step=10)
            #ode_runner.set_integrator('dop853', rtol=1e-6, atol=1e-6, max_step=10)
            ode_runner.set_initial_value(self.xis.copy()) 
            dt = self.par['t'] / self.par['Nt']
            i = 0
            self.z = zeros((self.Nt+2, self.Ntot), float)
            self.z[0, :] = self.xis[:]
            while ode_runner.successful() and ode_runner.t < self.t:
                yl = ode_runner.integrate(ode_runner.t + dt)
                i = i + 1
                self.z[i, :] = yl[:]
        else:
            self.z, self.infodict = odeint(self.funcT, self.xis, self.tv, full_output=True)
        sz = self.z.shape
        self.z = reshape(self.z, (sz[0], 2 + self.N_CMC, self.Np))
        self.z[:, :, 0] = 4. / 3. * self.z[:, :, 1] - 1. / 3. * self.z[:, :, 2]
        self.z[:, :, -1] = 4. / 3. * self.z[:, :, -2] - 1. / 3. * self.z[:, :, -3]
        #self.chi = self.itau * (self.z[:, 1, :] - self.z[:, 0, :]**2)


if __name__ == "__main__":
    #sol = mixinglayer(itau=10., L=8., t=0.2, Nt=100)
    sol = reactive_mixinglayer(vode=True, t=0.1, itau=10., L=4., N=40, Nt=200, A=5e7)

    sol.init()
    sol.solve()
    f  = sol.z[:, 0, :]
    im = sol.z[:, 1, :]
    
    s = im - f * f
    x = s / f / (1 - f)    
    contourf(x, linspace(0, 1, 20))
    colorbar()
    
    pmf = Conditionals(9, 25, sol, _plot=True)
    print pmf.fmean, pmf.sigma, pmf.chi, double_getitem(pmf.grad_f, 0), double_getitem(pmf.grad_s, 0)
    #nt, nz = x.shape
    #tau = zeros((nt, nz))
    #dtaudf = zeros((nt, nz))
    #dtauds = zeros((nt, nz))
    #pmf = None
    #for i in range(nt):
        #for j in range(nz):
            #_pmf = Conditionals(i, j, sol, pmfclass=pmf)
            #if _pmf:
                #pmf = _pmf
                #tau[i, j] = pmf.tau
                #dtaudf[i, j] = pmf.dtaudf
                #dtauds[i, j] = pmf.d2taudsds
    #show()
        
    #pmf = Conditionals(20, 20, sol, lookup="GSL_Lookuptable50.dat")
    #pmf = Conditionals(20, 20, sol)
    #show()
    