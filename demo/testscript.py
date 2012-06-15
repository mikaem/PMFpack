from PMFpack import *
from numpy import *
from pylab import *
from pygsl import integrate
from scipy.integrate import quad, dblquad, Inf
from scipy.special import erf, erfinv

w1 = integrate.workspace(10000)
w2 = integrate.workspace(10000)

#fmean = 0.639756248702123
#sigma = 0.16674922174125 * fmean * (1 - fmean)
#pmf = PMF(fmean, sigma, 0, 0, 0, 2)
#pmf.set_fmean_gradient(-0.4178939524532, 0, 0)  # Last two are zero because the simulation is 1D
#pmf.set_sigma_gradient(-0.5095137299633, 0, 0)
#pmf.chi = 0.38430391476813  # Mean scalar dissipation rate (defined with the factor 2)

#fmean = 0.32734785270559791
#sigma = 0.036258388156940768 * fmean * (1 - fmean)
#pmf = PMF(fmean, sigma, 0, 0, 0, 2)
#pmf.set_fmean_gradient(-0.40312611148071925, 0, 0)  # Last two are zero because the simulation is 1D
#pmf.set_sigma_gradient(-0.29357760186462173, 0, 0)
#pmf.chi = 0.36258388156940158  # Mean scalar dissipation rate (defined with the factor 2)

#fmean = 0.67265214729440204
#sigma = 0.036258388156940768 * fmean * (1 - fmean)
#pmf = PMF(fmean, sigma, 0, 0, 0, 2)
#pmf.set_fmean_gradient(-0.40312611148071925, 0, 0)  # Last two are zero because the simulation is 1D
#pmf.set_sigma_gradient(-0.51267462109681816, 0, 0)
#pmf.chi = 0.36258388156940158  # Mean scalar dissipation rate (defined with the factor 2)

fmean = 0.49994238
sigma = 0.13800082515647047
pmf = PMF(fmean, sigma, 0,0,0,2)
pmf.set_fmean_gradient(-2.22233706e+00,  -5.90368198e-05, 0)
pmf.set_sigma_gradient(-2.22227923, -0.30140048, 0)
pmf.chi = 50 * sigma

pmf.DT = 1.         # Turbulent diffusivity

pmf.compute(0, True)

def pdf(x, a, t):
    er = erfinv(2 * x - 1)
    return sqrt(2. * t / (1. - 2. * t)) * exp(er**2 - (a + 2 * sqrt(t) * er)**2 / 2 / (1 - 2 * t) )

def r(phi, a, t):
    return 1. / sqrt(2. * pi * (1. - 2. * t)) * exp(-phi**2 / 2. / (1. - 2. * t) )
    
def xpdf(x, a, t):
    er = erfinv(2*x-1)
    return x * pdf(x, a, t)

def xr(phi, a, t):
    xx = 0.5 * (1 + erf((phi - a) / 2. / sqrt(t)))
    return xx * r(phi, a, t)
        
def X(phi, a, t):
    return 0.5 * (1 + erf((phi - a) / 2. / sqrt(t)))

def phi(x, a, t):
    return a + 2 * sqrt(t) * erfinv(2 * x - 1)

def cdf(x, (a, t)):
    sys = integrate.gsl_function(pdf, (a, t))
    flag, result, error = integrate.qagp(sys, [0., x], 1e-8, 1e-7, 100, w1)
    return result
 
def cdf_eta(x, a, t):
    return quad(pdf, 1e-10, x, args=(a, t))[0]

def cdf_phi(x, a, t):
    return quad(r, -12., x, args=(a, t))[0]
    
def II(x, a, t, method=5):
    if method == 0:
        sys = integrate.gsl_function(cdf, (a, t))
        integral = integrate.qagp(sys, [0., x], 1e-8, 1e-7, 100, w2)
    elif method == 2:
        integral = quad(cdf_eta, 1e-8, x, args=(a, t))
    elif method == 3:
        integral = x * cdf_eta(x, a, t) - quad(xpdf, 1e-10, x, args=(a, t))[0] 
    elif method == 4:
        phi = a + 2 * sqrt(t) * erfinv(2. * x - 1.)
        integral = x * cdf_phi(phi, a, t) - quad(xr, -12., phi, args=(a, t))[0]
    elif method == 5:
        phi = a + 2 * sqrt(t) * erfinv(2. * x - 1.)
        integral = x / 2. * (1. + erf(phi / sqrt(2. * (1. - 2. * t)))) - quad(xr, -Inf, phi, args=(a, t))[0] 
    return integral

def dalfadf(a):
    return -sqrt(2 * pi) * exp(a * a / 2.)
    
def dpdfdtau(x, pmf):
    a, t = pmf.alfa, pmf.tau
    phi = a + 2 * sqrt(t) * erfinv(2 * x - 1)
    return pdf(x, a, t) / 2. / t / (1 - 2 * t) * (1 + phi * a - phi**2 / (1 - 2 * t))

def dpdfdtau_numeric(x, pmf, dh=1e-4):
    a, t = pmf.alfa, pmf.tau
    return (pdf(x, a, t+dh) - pdf(x, a, t-dh)) / 2 / dh 
    
def dpdfdalfa(x, pmf):
    a, t = pmf.alfa, pmf.tau
    return - pdf(x, a, t) * phi(x, a, t) / (1 - 2 * t)

def dpdfdalfa_numeric(x, pmf, dh=1e-4):
    a, t = pmf.alfa, pmf.tau
    return (pdf(x, a+dh, t) - pdf(x, a-dh, t)) / 2 / dh

def N(x, a, t):
    return sqrt((1 - t) / t) * exp(-2 * erfinv(2 * x - 1)**2 + a**2 / 2 / (1 - t))

def dNdtau(x, pmf):
    a, t = pmf.alfa, pmf.tau
    phi = a + 2 * sqrt(t) * erfinv(2 * x - 1)
    return N(x, a, t) * (-1. / (2. * t * (1 - t)) + a**2 / 2 / (1 - t)**2)

def dNdtau_numeric(x, pmf, dh=1e-4):
    a, t = pmf.alfa, pmf.tau
    return (N(x, a, t+dh) - N(x, a, t-dh)) / 2 / dh
    
def dNdalfa(x, pmf):
    a, t = pmf.alfa, pmf.tau
    return N(x, a, t) * a / (1 - t)
    
def dNdalfa_numeric(x, pmf, dh=1e-4):
    a, t = pmf.alfa, pmf.tau
    return (N(x, a+dh, t) - N(x, a-dh, t)) / 2 / dh
    
def dIIdtau(x, pmf):
    a, t = pmf.alfa, pmf.tau    
    return 0.5 * N(x, a, t) * pdf(x, a, t) / pmf.dtauds

def dIIdtau_numeric(x, pmf, dh=0.0001):
    a, t = pmf.alfa, pmf.tau
    return (II(x, a, t+dh) - II(x, a, t-dh)) / 2. / dh

def dIIdalfa(x, pmf):
    a, t = pmf.alfa, pmf.tau
    ex = erfinv(2 * x - 1)
    return - 1. / 2. / dalfadf(a) * (1 + erf((ex + sqrt(t) * a) / sqrt(1 - 2 * t)))

def dIIdalfa_numeric(x, pmf, dh=0.0001):
    a, t = pmf.alfa, pmf.tau
    return (II(x, a+dh, t) - II(x, a-dh, t)) / 2 / dh

def d2IIdalfa2(x, pmf):
    a, t = pmf.alfa, pmf.tau
    return exp(-a**2 / 2.) / 2. / sqrt(2.) * (2 * sqrt(t) / pi / sqrt(1 - 2 * t) * exp(-((erfinv(2 * x - 1) + a * sqrt(t)) / (sqrt(1 - 2 * t)))**2) - a / sqrt(pi) * (1 + erf((erfinv(2 * x - 1) + sqrt(t) * a) / (sqrt(1 - 2 * t)))))

def d2IIdalfa2_numeric(x, pmf, dh=0.0001):
    a, t = pmf.alfa, pmf.tau
    return (II(x, a+dh, t) - 2. * II(x, a, t) + II(x, a-dh, t)) / dh**2

def d2IIdtau2(x, pmf):
    a, t = pmf.alfa, pmf.tau    
    f1 = 0.5 / pmf.dtauds * (N(x, a, t) * dpdfdtau(x, pmf) + pdf(x, a, t) * dNdtau(x, pmf))
    f2 = 0.5 * N(x, a, t) * pdf(x, a, t) * pmf.d2taudsds / pmf.dtauds**3
    return f1 - f2

def d2IIdtau2_numeric(x, pmf, dh=0.0001):
    a, t = pmf.alfa, pmf.tau
    return (II(x, a, t+dh) - 2. * II(x, a, t) + II(x, a, t-dh)) / dh**2
    
def d2IIdadtau(x, pmf):
    a, t = pmf.alfa, pmf.tau
    phi = a + 2 * sqrt(t) * erfinv(2 * x - 1)
    return - phi / 2. / dalfadf(a) / sqrt(t * pi) / (1 - 2 * t)**1.5 * exp(-((erfinv(2 * x - 1) + sqrt(t) * a) / sqrt(1. - 2. * t))**2)
    
def d2IIdadtau_numeric(x, pmf, dh=0.0001):
    a, t = pmf.alfa, pmf.tau
    return (II(x, a+dh, t+dh) - II(x, a+dh, t-dh) - II(x, a-dh, t+dh) + II(x, a-dh, t-dh)) / 4. / dh**2
        
def dIIds(x, pmf):
    return 0.5 * N(x, pmf.alfa, pmf.tau) * pdf(x, pmf.alfa, pmf.tau)

def dIIds_numeric(x, pmf, dh=1e-5):
    pmf.set_parameters(pmf.fmean, pmf.sigma + dh)
    pmf.compute(0, False)
    ip = II(x, pmf.alfa, pmf.tau)
    pmf.set_parameters(pmf.fmean, pmf.sigma - 2 * dh)
    pmf.compute(0, False)
    im = II(x, pmf.alfa, pmf.tau)
    pmf.set_parameters(pmf.fmean, pmf.sigma + dh)
    pmf.compute(0, False)
    return (ip - im) / 2 / dh
    
def d2IIdsds(x, pmf):
    f1 = N(x, pmf.alfa, pmf.tau) * dpdfdtau(x, pmf)
    f2 = pdf(x, pmf.alfa, pmf.tau) * dNdtau(x, pmf)
    return 0.5 * (f1 + f2) * pmf.dtauds

def d2IIdsds_numeric(x, pmf, dh=1e-5):
    pmf.set_parameters(pmf.fmean, pmf.sigma + dh)
    pmf.compute(0, False)
    ip = II(x, pmf.alfa, pmf.tau)
    pmf.set_parameters(pmf.fmean, pmf.sigma - 2 * dh)
    pmf.compute(0, False)
    im = II(x, pmf.alfa, pmf.tau)
    pmf.set_parameters(pmf.fmean, pmf.sigma + dh)
    pmf.compute(0, False)
    ii = II(x, pmf.alfa, pmf.tau)
    return (ip - 2. * ii + im) / dh**2
    
def dIIdf(x, pmf, central=False):
    pmf.central = central
    pmf.compute(0, True)
    return dIIdalfa(x, pmf) * dalfadf(pmf.alfa) + dIIdtau(x, pmf) * pmf.dtaudf

def dIIdf_numeric(x, pmf, dh=1e-5, central=False):
    f0, s0, im0 = pmf.fmean, pmf.sigma, pmf.im
    if not central:
        ds = - 2. * f0 * dh - dh**2
        pmf.set_parameters(f0 + dh, s0 + ds)
        pmf.compute(0, False)
        ip = II(x, pmf.alfa, pmf.tau)
        ds = 2. * f0 * dh - dh**2
        pmf.set_parameters(f0 - dh, s0 + ds)
        pmf.compute(0, False)
        im = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0, s0)
        pmf.compute(0, False)
    else:
        pmf.set_parameters(f0 + dh, s0)
        pmf.compute(0, False)
        ip = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0 - dh, s0)
        pmf.compute(0, False)
        im = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0, s0)
        pmf.compute(0, False)
        
    return (ip - im) / 2 / dh
    
def d2IIdfds(x, pmf, central=False):
    pmf.central = central
    pmf.compute(0, True)
    a, t = pmf.alfa, pmf.tau
    f1 = N(x, a, t) * (dpdfdalfa(x, pmf) * dalfadf(a) + dpdfdtau(x, pmf) * pmf.dtaudf)
    f2 = pdf(x, a, t) * (dNdalfa(x, pmf) * dalfadf(a) + dNdtau(x, pmf) * pmf.dtaudf)
    return 0.5 * (f1 + f2)
    
def d2IIdfds_numeric(x, pmf, dh=1e-5, central=False):
    pmf.central = central
    f0, s0 = pmf.fmean, pmf.sigma
    if not central:
        ds = (1. - 2. * f0) * dh - dh**2
        pmf.set_parameters(f0 + dh, s0 + ds)
        pmf.compute(0, False)
        ipp = II(x, pmf.alfa, pmf.tau)
        ds = (1. + 2. * f0) * dh - dh**2
        pmf.set_parameters(f0 -  dh, s0 + ds)
        pmf.compute(0, False)
        ipm = II(x, pmf.alfa, pmf.tau)
        ds = -(1. + 2. * f0) * dh - dh**2
        pmf.set_parameters(f0 + dh, s0 + ds)
        pmf.compute(0, False)
        imp = II(x, pmf.alfa, pmf.tau)
        ds = (-1. + 2. * f0) * dh - dh**2
        pmf.set_parameters(f0 - dh, s0 + ds)
        pmf.compute(0, False)
        imm = II(x, pmf.alfa, pmf.tau)
    else:
        pmf.set_parameters(f0 + dh, s0 + dh)
        pmf.compute(0, False)
        ipp = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0 -  dh, s0 + dh)
        pmf.compute(0, False)
        ipm = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0 + dh, s0 - dh)
        pmf.compute(0, False)
        imp = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0 - dh, s0 - dh)
        pmf.compute(0, False)
        imm = II(x, pmf.alfa, pmf.tau)
    pmf.set_parameters(f0, s0)
    pmf.compute(0, False)        
    return (ipp + imm - ipm - imp) / 4 / dh / dh

def d2IIdfdf(x, pmf, central=False):
    pmf.central = central
    pmf.compute(0, True)
    a, t = pmf.alfa, pmf.tau
    f1 = - dalfadf(a) * sqrt(t / pi / (1. - 2. * t)) * exp(-((erfinv(2 * x - 1) + sqrt(t) * a) / sqrt(1. - 2. * t))**2)
    #f1 = dalfadf(a)**2 * (d2IIdalfa(x, a, t) + a * dIIdalfa(x, a, t))
    f2 = 2. * d2IIdadtau(x, pmf) * dalfadf(a) * pmf.dtaudf
    f3 = d2IIdtau2(x, pmf) * pmf.dtaudf**2
    f4 = dIIdtau(x, pmf) * pmf.d2taudfdf
    return f1 + f2 + f3 + f4

def d2IIdfdf_numeric(x, pmf, dh=1e-5, central=False):
    f0, s0, im0 = pmf.fmean, pmf.sigma, pmf.im
    if not central:
        ds = - 2. * f0 * dh - dh**2
        pmf.set_parameters(f0 + dh, s0 + ds)
        pmf.compute(0, False)
        ip = II(x, pmf.alfa, pmf.tau)
        ds = 2. * f0 * dh - dh**2
        pmf.set_parameters(f0 - dh, s0 + ds)
        pmf.compute(0, False)
        im = II(x, pmf.alfa, pmf.tau)
    else:
        pmf.set_parameters(f0 + dh, s0)
        pmf.compute(0, False)
        ip = II(x, pmf.alfa, pmf.tau)
        pmf.set_parameters(f0 - dh, s0)
        pmf.compute(0, False)
        im = II(x, pmf.alfa, pmf.tau)        
    pmf.set_parameters(f0, s0)
    pmf.compute(0, False)
    ii = II(x, pmf.alfa, pmf.tau)
    return (ip - 2. * ii + im) / dh**2

def CSD(x, pmf):
    df1 = double_getitem(pmf.grad_f, 0)
    ds1 = double_getitem(pmf.grad_s, 0)
    df2 = double_getitem(pmf.grad_f, 1)
    ds2 = double_getitem(pmf.grad_s, 1)
    fx2 = df1 * df1 + df2 * df2
    sx2 = ds1 * ds1 + ds2 * ds2
    fxsx = df1 * ds1 + df2 * ds2
    a, t = pmf.alfa, pmf.tau
    f1 = N(x, a, t) * pmf.chi
    f2 = 2. * pmf.DT / pdf(x, a, t) * (
        d2IIdsds(x, pmf) * sx2 +
        d2IIdfdf(x, pmf) * fx2 +
        2. * d2IIdfds(x, pmf) * fxsx
    )
    return f1 + f2

def CSD_numeric(x, pmf):
    df1 = double_getitem(pmf.grad_f, 0)
    ds1 = double_getitem(pmf.grad_s, 0)
    df2 = double_getitem(pmf.grad_f, 1)
    ds2 = double_getitem(pmf.grad_s, 1)
    fx2 = df1 * df1 + df2 * df2
    sx2 = ds1 * ds1 + ds2 * ds2
    fxsx = df1 * ds1 + df2 * ds2
    a, t = pmf.alfa, pmf.tau
    f1 = N(x, a, t) * pmf.chi
    f2 = 2. * pmf.DT / pdf(x, a, t) * (
        d2IIdsds_numeric(x, pmf, dh=1e-5) * sx2 +
        d2IIdfdf_numeric(x, pmf, dh=1e-5) * fx2 +
        2. * d2IIdfds_numeric(x, pmf, dh=1e-5) * fxsx
    )
    return f1 + f2

tests = ['dpdfdtau', 'dpdfdalfa', 
         'dNdtau', 'dNdalfa', 
         'dIIdalfa', 'dIIdtau',
         'd2IIdtau2', 'd2IIdalfa2', 'd2IIdadtau',
         'dIIds', 'dIIdf', 'd2IIdsds', 'd2IIdfdf', 'd2IIdfds',
         'CSD']
x = 0.5
for test in tests:
    exec('a = ' + test + '(x, pmf); b = ' + test + '_numeric(x, pmf)')
    print test, ' error = ', a - b, ' vals ', a, b
print pmf.CSD(x), pmf.CSD_verify(x), pmf.counterflow(x)
