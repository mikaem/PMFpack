"""
We study non-premixed combustion in a turbulent mixing layer using a 
simplified reaction model

    F + r * O  <--> (1 + r) * P
    
where non-premixed Fuel (F) and Oxidizer (O) combust to form Product 
(P). r is the stoichiometric ratio or the mass of oxidant disappearing 
with one unit mass of fuel.

The mixture fraction, xi, is a conserved scalar that can be related
to the mass fractions of the reacting fuel (YF) and oxidizer (YO) 
through

    xi = (r * YF - YO + 1) / (r + 1)
   
This simplified combustion problem can be solved using only the 
mixture fraction and one additional parameter, the reduced
temperature, T, related to fuel and oxidizer through

    YF = xi - T / (r + 1)
    YO = 1 - xi - r * T / (r + 1)

We solve an equation for the conserved mixture fraction mean (<xi>) 
and second raw moment (<xi*xi>), given by the equations

    U grad(<xi>)    = div(grad(<xi>))    
    U grad(<xi*xi>) = div(grad(<xi*xi>)) - itau * ( <xi*xi> - <xi>^2 )

where U here is set to the constant vector (0, 10).

The square domain is x = [-L/2, L/2], y = [0, 8] with Neumann boundary 
conditions at the ends. The inlet conditions at y = 0 are
    
    <xi> = 1 - H(x), <xi*xi> = <xi>,
    
where H(x) is the Heaviside function that is unity for 
positive x and zero otherwise. For stability the Heviside function is 
approximated with an error function.
    
The central variance is given as sigma = <xi*xi> - <xi>^2

The parameter itau determines the rate of mixing compared to transport.
The mean scalar dissipation rate (chi) is given as

    chi = itau * (<xi*xi> - <xi>^2)
    
The instantaneous reaction rate is given as

    rx = (r + 1) A exp(- Ze/alpha - Ze (1 - T) / (1 - alpha (1 - T))) * 
         (YF * YO^r - kappa * T^(r + 1))

where Ze, A, alpha and kappa are parameters and T is the reduced 
temperature.

We use three different models for the combustion, where the most
advanced is the conditional moment closure (CMC). With CMC we solve for 
the mean temperature conditional on the value of xi, i.e., <T|eta>, 
where eta is a sample space variable for xi. The CMC equation for 
<T|eta> is (setting T = <T|eta> for simplicity of notation)

    (U + 2 <u'|eta>) grad(T) = div(grad(T)) + 0.5 <chi|eta> d^2T / deta^2 + rx(eta, T)

where <u'|eta> and <chi|eta> are the conditional mean velocity and
conditinal scalar dissipation rate respectively.

The steady flamelet model is retained by using only the last two terms 
on the right hand side, i.e., no explicit spatial dependency for <T|eta>.

The equilibrium model is retained by setting the reaction rate equal
to zero.

The unconditional temperature can be obtained from

    <T> = \int_0^1 pdf(eta) * <T|eta> deta

and the unconditional concentrations of fuel and oxidizer are found
similarily using the above expressions for YF and YO.

The spatial dependency of <T|eta> is weak, so we can solve the CMC
equations on a coarser mesh than the unconditional <xi> and <xi*xi>.
We do this here by using a high order approximation for the unconditional
moments <xi> and <xi*xi>, and a first order approximation for <T|eta>.

The CMC equation is solved using finite elements in space and a central 
finite difference stencil in mixture fraction space. We solve the 
equation using full coupling in mixture fraction space and a fully non-
linear iterative Newton method.

Conditional models for <u'|eta>, <chi|eta> and the pdf are closed using 
the presume mapping function (PMF) approach, provided by the PMFpack 
module.
"""
from PMFpack import *
from dolfin import *
#from cbc.pdesys import *
from numpy import array, isnan, sqrt, linspace, sign, cos, zeros, maximum, minimum
from pylab import find
from scipy.special import erf
from scipy.optimize import fsolve
from time import time

# Solve unconditional models for flim < fmean < 1 - flim
# Assume equilibrium otherwise
flim = 0.005

# Create a global pmf class for computing conditional models
pmf = PMF(0.5, 0.1, 0, 0, 0, 0)

# Reaction rate constants
r = Constant(3)
kappa = Constant(0.01)
A = Constant(1.0e7)
alfa = Constant(0.87)
Ze = Constant(4.)

eq_solution = {}
class ReactionModel:
    """Compute the steady flamelet or equilibrium model for given parameters.
    This is computed in seperate class since for these models there is no 
    explicit spatial dependency for <T|eta>."""
    def __init__(self, eta):
        self.eta = eta
        self.mesh = Interval(100, 0, 1)
        V = FunctionSpace(self.mesh, 'CG', 1)
        u = TrialFunction(V)
        v = TestFunction(V)
        self.u_ = u_ = Function(V)
        self.csd = csd = Function(V)
        self.csd_array = zeros(V.dim())
        
        # Define flamelet problem
        eta = Expression("x[0]")
        form = inner(grad(u_), grad(dot(v, csd / 2.)))*dx - inner(v, 
            (r + 1) * A * exp(- Ze / alfa - Ze * (1 - u_) / (1 - alfa * (1 - u_))) * 
            ((eta - u_ / (r + 1)) * pow((1 - eta - r * u_ / (r + 1)), r) - kappa * pow(u_, r + 1)))*dx
        
        # Create boundary conditions
        self.bc = DirichletBC(V, Constant(0), "x[0] < DOLFIN_EPS || x[0] > 1 - DOLFIN_EPS")
        
        # Create Jacobian and declare solution tensors
        self.J = derivative(form, u_, u)
        self.F = form
        self.y = interpolate(eta, V)
        self.yy = yy = self.y.vector().array()
        self.eta_valid_index = find(yy * (1. - yy) > 0.)
        self.y_valid = self.yy[self.eta_valid_index]
        self.cx = zeros(len(self.eta_valid_index))
        
    def flamelet(self, fmean, sigma, chi):
        """Solve steady flamelet model using Newton iterations."""
        iseg = sigma / fmean / (1 - fmean)
        self.EQ() # Initialize Newton solver, or use this solution for extremes
        if (fmean < flim or fmean > (1 - flim) or  iseg < flim or iseg > (1 - flim)):
            pass
        else:
            pmf.set_parameters(fmean, sigma)
            pmf.chi = chi
            pmf.compute(0, False)
            pmf.counterflow(self.y_valid, self.cx)
            self.csd_array[self.eta_valid_index] = self.cx[:]
            self.csd.vector().set_local(self.csd_array)
            solve(self.F == 0, self.u_, [self.bc], J=self.J, solver_parameters={"newton_solver":{'report': False}})
            
    def EQ(self):
        """Equilibrium solution."""
        def rx(x, eta):
            """Reaction rate."""
            return (eta - x / (r(0)+1)) * (1 - eta - r(0) * x / (r(0) + 1))**r(0) - kappa(0) * x**(r(0)+1)

        for i, eta in enumerate(self.yy):
            if i == 0 or i == len(self.yy)-1:
                self.u_.vector()[i] = 0.
            elif r(0) == 1: # Second degree polynomial
                a = r(0) / (r(0) + 1.)**2 - kappa(0)
                b = - 1. / (r(0) + 1.)
                c = eta - eta**2
                self.u_.vector()[i] = (-b - sqrt(b**2 - 4 * a * c)) / 2. / a
            else:
                if (eta, r(0), kappa(0)) in eq_solution:
                    self.u_.vector()[i] = eq_solution[(eta, r(0), kappa(0))]
                else:
                    val = fsolve(rx, eta*(1-eta), args=(eta))
                    self.u_.vector()[i] = val[0]
                    eq_solution[(eta, r(0), kappa(0))] = val[0]

class Tinit(Expression):
    """Compute equilibrium or flamelet solution."""    
    def __init__(self, eta, fm=None, model=0, itau=1., **kwargs):
        self.rx = ReactionModel(eta)
        self.model = model
        self.fm = fm
        self.itau = itau
        self.computed_sol = {} # This is used to hold the solution in the nodes to avoid costly recomputations
        
    def eval(self, values, x):
        tx = tuple(x)
        if tx in self.computed_sol:
            sol = self.computed_sol[tx]
        else:    
            if self.model == 0:
                self.rx.EQ()                
            elif self.model == 1:
                f, im = self.fm(x)
                s = im - f * f
                chi = s * self.itau(x)
                self.rx.flamelet(f, s, chi)
            sol = array([self.rx.u_(xx) for xx in self.rx.eta[1:-1]])
            self.computed_sol[tx] = sol
            
        values[:] = sol[:]
                    
class Conditionals(Expression):
    """Compute CSD and conditional velocity in x- and y-direction.
    """
    def __init__(self, Nc, fm=None, grad_f=None, grad_im=None, itau=1, **kwargs):
        self.fm = fm             # mean and second integer moment
        self.grad_f = grad_f
        self.grad_im = grad_im
        self.itau = itau
        self.csd = zeros(Nc)     # Conditional scalar dissipation rate
        self.cv = zeros((3, Nc)) # Conditional velocity
        self.computed_sol = {}
    
    def eval(self, values, x):
        tx = tuple(x)
        N = len(self.csd)
        self.csd[:] = 0
        self.cv[:] = 0
        if tx in self.computed_sol:
            self.csd[:], self.cv[:] = self.computed_sol[tx]
        else:
            f, im = self.fm(x)
            f = max(min(1., f), 0.)
            s = max(0., min(im - f * f, f * (1. - f)))
            chi = self.itau(x) * s
            if f > flim and f < (1 - flim):
                seg = s / f / (1 - f)
                if seg > flim and seg < (1 - flim):
                    pmf.set_parameters(f, s)
                    pmf.set_fmean_gradient(*self.grad_f(x))
                    pmf.set_sigma_gradient(*self.grad_im(x))
                    pmf.chi = chi
                    pmf.compute(0, True)
                    pmf.CV(eta[1:-1], self.cv)
                    if seg > 0.95:
                        pmf.counterflow(eta[1:-1], self.csd)
                    else:
                        pmf.counterflow(eta[1:-1], self.csd)
                        #pmf.CSD(eta[1:-1], self.csd)
            self.computed_sol[tx] = self.csd.copy(), self.cv.copy()
            
        values[:N] = self.csd[:]
        values[N:2*N] = self.cv[0, :]
        values[2*N:3*N] = self.cv[1, :]
                
class Finit(Expression):
    """Initialize mean mixture fraction and second integer moment."""
    def __init__(self, sigma=0.01, **kwargs):
        self.sigma = sigma
        
    def eval(self, values, x):
        values[0] = (1 - 0.5 * (1 + erf(x[0] / self.sigma)))
        values[1] = values[0]
        
    def value_shape(self):
        return (2,)
        
def mixinglayermesh(L=2, H=8, Nx=40, Ny=40):
    mesh = Rectangle(-L, 0, L, H, Nx, Ny)
    x = mesh.coordinates()
    x[:, 0] = sign(x[:, 0]) * (1. - cos(pi * x[:, 0] / L / 2.)) * L
    return mesh
    
# Set up problem with mesh and function spaces
mesh = mixinglayermesh(Nx=40, Ny=40)
V = FunctionSpace(mesh, 'CG', 1) 
FM = MixedFunctionSpace([V, V])  # fmean, variance
fm = TrialFunction(FM)
fm_v = TestFunction(FM)
f0 = interpolate(Finit(), FM)
fm_ = f0.copy(True)
x = Vector(f0.vector())

# Dirichlet boundary condition on inlet
bc = DirichletBC(FM, f0, "on_boundary && x[1] < 10 * DOLFIN_EPS")

# Normalized velocity
U0 = Constant((0., 10.))

# Normalized mixing frequency (Need this to be zero at x[1] = 0 to be consistent with eq-solution. If flamelet solution is used as Dirichlet on inlet it can be nonzero.)
#itau = Constant(10.)
itau_e = Expression("10. * (1. - exp(-5. * x[1] * x[1] * x[1]))")
itau = interpolate(itau_e, V)

# Variational form for first two integer moments of mixture fraction
F0 = inner(grad(fm_) * U0, fm_v)*dx + inner(grad(fm_), grad(fm_v))*dx + itau * inner(fm_[1] - fm_[0] * fm_[0], fm_v[1])*dx

# Compute derivative of nonlinear form
J0 = derivative(F0, fm_, fm)

# Solve using Newton iterations
solve(F0 == 0, fm_, [bc], J=J0)
#scalar = PDESubSystem(vars(), ['fm'], bcs=[bc], F=F0, iteration_type='Newton')
#scalar.solve(max_iter=5)

# CMC equations are solved on top of the scalar moments
# Divide eta space into Nc-1 intervals
Nc = 10
CMCmesh = mixinglayermesh(Nx=10, Ny=10)
eta = linspace(0, 1, Nc)
dd = Constant((1. / (Nc - 1))**2)
Q = FunctionSpace(CMCmesh, 'CG', 1)
VQ = MixedFunctionSpace([Q,] * (Nc-2)) # Q(eta=1/Nc), Q(eta=2/Nc), ..., Q(eta=1-1/Nc)
#DQ = FunctionSpace(mesh, 'DG', 0)      # Alternatively one could use DQ for conditinal models
#DVQ = MixedFunctionSpace([DQ,] * (Nc-2))
q = TrialFunction(VQ)
q_v = TestFunction(VQ)

# Equilibrium solution
qeq = interpolate(Tinit(eta, fm=fm_, model=0, itau=itau, element=VQ.ufl_element()), VQ)

# Flamelet solution
qfl = interpolate(Tinit(eta, fm=fm_, model=1, itau=itau, element=VQ.ufl_element()), VQ)

# CMC solution Function q_
q_ = qfl.copy(True)

# Dirichlet boundary condition on inlet
bcmc = DirichletBC(VQ, qeq, "on_boundary && x[1] < 10 * DOLFIN_EPS")

# Compute conditional models
grad_f = project(grad(fm_[0]), Q * Q)
grad_im = project(grad(fm_[1]), Q * Q)
VVQ = MixedFunctionSpace([VQ, VQ, VQ]) # CSD + cond. velocity x and y direction
conditionals = Conditionals(len(eta)-2, fm=fm_, grad_f=grad_f, grad_im=grad_im, itau=itau, element=VVQ.ufl_element())
cmixed = interpolate(conditionals, VVQ)
CSD, cvx, cvy = cmixed.split()

print "Computed conditionals "

# Set up CMC-form for Newton iterations
def reaction(i):
    return - (inner( (r + 1) * A * exp(- Ze / alfa - Ze * (1 - q_[i]) / (1 - alfa * (1 - q_[i])))
           * ((eta[i+1] - q_[i] / (r + 1)) * pow((1 - eta[i+1] - r * q_[i] / (r + 1)), r(0)) - kappa * pow(q_[i], r(0) + 1)), q_v[i]) * dx)
           
CMC  =  inner(grad(q_) * U0, q_v)*dx + inner(grad(q_), grad(q_v))*dx

# Mixture fraction space uses finite difference discretization
for i in range(Nc-2):
    # Boundary conditions at i = 0 and i = Nc-3
    if i == 0:
        CMC -= CSD[0] / 2. * inner((q_[1] - 2 * q_[0]) / dd , q_v[0])*dx
    elif i == (Nc-3):
        CMC -= CSD[Nc-3] / 2. * inner((- 2 * q_[Nc-3] + q_[Nc-4]) / dd , q_v[Nc-3])*dx    
    else:
        CMC -= CSD[i] / 2. * inner((q_[i+1] - 2 * q_[i] + q_[i-1]) / dd , q_v[i])*dx 
    CMC += 2. * cvx[i] * inner(grad(q_[i])[0], q_v[i])*dx
    CMC += 2. * cvy[i] * inner(grad(q_[i])[1], q_v[i])*dx
    CMC += reaction(i)

# Compute derivative of nonlinear form
J1 = derivative(CMC, q_, q)   

# Solve CMC equation using Newton iterations
solve(CMC == 0, q_, [bcmc], J=J1)
#cmcsolver = PDESubSystem(vars(), ['q'], bcs=[bcmc], F=CMC, iteration_type='Newton')
#cmcsolver.solve(max_iter=5)

# Compute unconditional average temperature by summing pdf(eta) * Q(eta)
# This should be integrated much more accurately using the exact pdf and
# (linear) interpolation of the conditional temperature
VV = MixedFunctionSpace([V,] * (Nc - 2))
pdf = Function(VV)
fa = fm_.vector().array()
N = V.dim()
fa[:N] = maximum(minimum(fa[:N], 1.), 0.)
sa = maximum(0, minimum(fa[N:] - fa[:N] * fa[:N], fa[:N] * (1. - fa[:N])))
T = Function(V)   # CMC unconditional temperature
Tfl = Function(V) # Flamelet unconditional temperature
qq = interpolate(q_, VV)   # Interpolate CMC solution on mixture fraction mesh
qf = interpolate(qfl, VV)  # Interpolate flamelet solution on mixture fraction mesh
pdf_ = zeros(Nc - 2)
for i in range(V.dim()):
    f = fa[i]
    if f > flim and f < (1 - flim):
        s = sa[i]
        seg = s / f / (1 - f)
        if seg > flim and seg < (1 - flim):
            pmf.set_parameters(f, s)
            pmf.compute(0, False)
            pmf.PDF(eta[1:-1], pdf_)
            for j in range(Nc-2):
                pdf.vector()[j * V.dim() + i] = pdf_[j]
                T  .vector()[i] += pdf_[j] * qq.vector()[j * V.dim() + i] * eta[1]
                Tfl.vector()[i] += pdf_[j] * qf.vector()[j * V.dim() + i] * eta[1]

plot(T, title="CMC Temperature")
plot(Tfl, title="Flamelet Temperature")
