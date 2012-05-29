// Copyright (C) 2012 Mikael Mortensen
//
// This file is part of PMFpack.
//
// PMFpack is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PMFpack is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with PMFpack. If not, see <http://www.gnu.org/licenses/>.
//
//
// First added:  2012-05-04
// Last changed: 2012-05-04
#ifdef HAS_BOOST
#include <boost/math/special_functions/erf.hpp>
#endif
#include "PMF.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_cdf.h"

// The implemented error function and its inverse are defined using 
// the normal distribution and not the error function as defined in
// e.g. http://mathworld.wolfram.com/Erf.html. The correlation between 
// the regular error function and the implemented based on the normal 
// distribution is as folows
//
//  erf(z) = 0.5 * (1 + erf_(z/sqrt(2)))
//
//  erfinv(z) = sqrt(2) * erfinv_(2 * z - 1)
// 
//  where
//                              / x
//                     2       |       -t*t
//       erf_(x) = ----------- |      e       dt
//                 sqrt(pi)    |
//                             /0
//
//                              / x
//                     1       |       -t*t/2
//        erf(x) = ----------- |      e       dt
//                  sqrt(2 pi) |
//                             /-oo
// 
// and erfinv and erfinv_ are the inverses of erf and erf_ respectively

#ifdef HAS_BOOST
// Boost use the erf_ definition and as such we need to wrap these up
double boost_erf(double z)
{
  return 0.5 * (1 + boost::math::erf(z / M_SQRT2));
}

double boost_erfinv(double z)
{
  return M_SQRT2 * boost::math::erf_inv(2 * z - 1);
}
#endif
namespace pmfpack
{  
  // Choose factory for error functions. I usually find GSL to be fastest (M. M.)
  simple_function erfinv = &gsl_cdf_ugaussian_Pinv;
//  simple_function erfinv = &boost_erfinv;
//   simple_function erfinv = &normsinv;
//   simple_function erfinv = &Erfinv;
//   simple_function erf = &Erf;  
  simple_function erf = &gsl_cdf_ugaussian_P;
//  simple_function erf = &boost_erf;
  
  void PMF::reallocate_solver(Root *root, int solver)
  {
    root->realloc(solver);
  }
  
  PMF::PMF(double _fmean, double _sigma, int fsolver, int fdfsolver, int integral, int deriv)
  {
    central = false; // Are you using central moments (or integer) in your CFD calculations?
    roots = new Roots();
    roots->froot = NULL;    
    roots->fdfroot = NULL;    
    roots->central = &central; 
    integrator = NULL;    
    derivator = NULL;
    verbose = 0;    
    chi = 1.;    
    DT = 1.;
    cheb_order = 6;
    set_parameters(_fmean, _sigma);
    parameters[0] = &fmean;
    parameters[1] = &sigma;
    parameters[2] = &alfa;
    parameters[3] = &tau;
    parameters[4] = &im;
    parameters[5] = &f;
    parameters[6] = &df;
    parameters[7] = &dtaudf;
    parameters[8] = &dtauds;
    parameters[9] = &d2taudfdf;
    parameters[10] = &d2taudsds;
    
    // Hook up the appropriate classe to do calculations
    set_integrator(integral);         
    set_fsolver(fsolver);     
    set_fdfsolver(fdfsolver);    
    set_derivator(deriv);
  }
  
  PMF::~PMF()
  {
    delete roots->froot;
    delete roots->fdfroot;
    delete roots;
    delete integrator;
    delete derivator;
  }
  
  void PMF::set_fsolver(int fsolver)
  {
    if (roots->froot != NULL) delete roots->froot;    
    
    if (fsolver == 0)
    {
      std::cout << "GSL fsolver" << std::endl; 
      roots->froot = new GSL_F_Root(integrator);
    }
#ifdef HAS_BOOST    
    else if (fsolver == 1)
    {
      std::cout << "Boost fsolver" << std::endl; 
      roots->froot = new Boost_F_Root(integrator);
    }
#endif
    else
    {
      std::cout << " Not implemented fsolver \n" << std::endl;
    }
  }
  
  void PMF::set_fdfsolver(int fdfsolver)
  {
    if (roots->fdfroot != NULL) delete roots->fdfroot;    
    
    if (fdfsolver == 0)
    {
      std::cout << "GSL fdfsolver" << std::endl;
      roots->fdfroot = new GSL_FDF_Root(integrator);
    }
#ifdef HAS_BOOST
    else if (fdfsolver == 1)
    {
      std::cout << "Boost fdfsolver" << std::endl; 
      roots->fdfroot = new Boost_FDF_Root(integrator);
    }
#endif
    else
    {
      std::cout << " Not implemented fdfsolver \n" << std::endl;
    }
  }
  
  void PMF::set_integrator(int integral)
  {
    if (integrator != NULL) delete integrator;
    
    if (integral == 0)
    {
      std::cout << "GSL integrator" << std::endl; 
      integrator = new GSLIntegrator(parameters);
    }
    else
    {
      std::cout << " Not implemented integrator \n" << std::endl;
    }
  }
  
  void PMF::set_derivator(int deriv)
  {
    if (derivator != NULL) delete derivator;
    
    if (deriv == 0)
    {
      std::cout << "Finite difference derivator" << std::endl; 
      derivator = new FDDerivator(roots);
    }
    else if (deriv == 1)
    {
      std::cout << "Chebyshev derivator" << std::endl; 
      derivator = new ChebDerivator(roots, cheb_order);
    }
    else
    {
      std::cout << " Not implemented derivator \n" << std::endl;
    }
  }
  
  void PMF::set_parameters(double f, double s)
  {
    fmean = f;
    sigma = s;
    alfa = erfinv(1 - f);
    tau = 0.1;
    im = s + f * f;
    if(f <= 0. || f >= 1.)
    {
      std::cout << "Wrong input. 0 < fmean < 1." << std::endl;
    }
    else if(s <= 0. || s / f / (1 - f) >= 1.){
      std::cout << "Wrong input. 0 < Intensity of segregation < 1." << std::endl;
    }
  }
  
  void PMF::set_fmean_gradient(double dfdx, double dfdy, double dfdz)
  {
    grad_f[0] = dfdx;
    grad_f[1] = dfdy;
    grad_f[2] = dfdz;
  }

  void PMF::set_sigma_gradient(double dsdx, double dsdy, double dsdz)
  {
    grad_s[0] = dsdx;
    grad_s[1] = dsdy;
    grad_s[2] = dsdz;
  }
  
  void PMF::compute(int solver, bool derivatives)
  {
    if (!derivatives){ // Compute just tau
      if (solver == 0)
        roots->froot->compute(verbose);     // Bracketing solver
      else if (solver == 1)
        roots->fdfroot->compute(verbose);   // Root polishing solver
      else
        std::cout << " Not implemented PMF::compute " << std::endl;
    }
    else{ // Compute tau, dtaudf, dtauds, d2taudfdf and d2taudsds
      derivator->compute(verbose);      
    }
  }
  
  double PMF::PDF(double eta)
  {
    double SQRT2tau, phi, EXV;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    SQRT2tau = sqrt(2. * tau);
    phi = alfa + SQRT2tau * EXV;
    return  SQRT2tau / sqrt(1. - 2. * tau) * exp(EXV * EXV / 2. - phi * phi / (2. * (1. - 2. * tau)));
  }

  void PMF::PDF(double *eta, int N1, double *p, int N2)
  {
    double SQRT2tau;
    double *phi, *EXV;

    assert(N1 == N2);
    EXV = (double *) malloc(N1 * sizeof(double));
    phi = (double *) malloc(N1 * sizeof(double));
    SQRT2tau = sqrt(2. * tau);
    for (int i = 0; i < N1 ; i++){
      EXV[i] = erfinv(eta[i]);
      phi[i] = alfa + SQRT2tau * EXV[i];
      p[i] = SQRT2tau / sqrt(1. - 2. * tau) * exp(EXV[i] * EXV[i] / 2. - phi[i] * phi[i] / (2. * (1. - 2. * tau)));
    }
    free(EXV);
    free(phi);
  }
  
  double PMF::counterflow(double eta)
  {
    double tm, EXV;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    tm = 1. - tau;
    return chi * sqrt(tm / tau) * exp(- EXV * EXV + alfa * alfa / 2. / tm);
  }
  
  void PMF::counterflow(double *eta, int N1, double *p, int N2)
  {
    double tm;
    double *EXV;

    assert(N1 == N2);
    EXV = (double *) malloc(N1 * sizeof(double));
    tm = 1. - tau;
    for (int i = 0; i < N1 ; i++){
      EXV[i] = erfinv(eta[i]);
    }
    for (int i = 0; i < N1 ; i++){
      p[i] = chi * sqrt(tm / tau) * exp(-EXV[i] * EXV[i] + alfa * alfa / 2. / tm);
    }
    free(EXV);
  }
  
  double PMF::CSD(double eta)
  {
    double alfaf, AA, BB, CC, CF, fx2, sx2, fxsx;
    double R1, R2, R3, R4, R5, SQRTtm2, SQRTtm, SQRTtau, SQRT2t, SQRT2tau, ALFA2;
    double tm, tm2, P1, cen;    
    double EXV, EXV2, phi, wa, csd;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    cen = central ? 2. : 0.;
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    SQRTtau = sqrt(tau);
    SQRT2t = SQRTtau * M_SQRT2;
    SQRTtm = sqrt(tm);
    SQRTtm2 = sqrt(tm2);
    R1 = SQRTtau / SQRTtm;
    R2 = SQRT2t / SQRTtm2;
    R3 = SQRTtm / SQRTtm2;
    SQRT2tau = SQRT2t * SQRTtm2;
    ALFA2 = alfa * alfa;
    alfaf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    fx2  = cblas_ddot(3, grad_f, 1, grad_f, 1);
    sx2  = cblas_ddot(3, grad_s, 1, grad_s, 1);
    fxsx = cblas_ddot(3, grad_f, 1, grad_s, 1);
    EXV2 = EXV * EXV;
    phi = alfa + SQRT2t * EXV;
    wa = exp(-EXV2);
    R4 = exp(ALFA2 / 2.);
    R5 = exp(ALFA2 / 2. / tm);
    BB = ALFA2 / 2. / tm /tm - (phi * EXV / M_SQRT2 / SQRTtau + phi * phi / tm2 - 0.5 / tm) / tm2;
    CF = SQRTtm / SQRTtau * wa * R5;
    AA = alfa / tm - phi / tm2;
    CC = M_SQRT2 / M_SQRTPI * (alfaf + phi * dtaudf / tau / tm2) * wa * R4;
    csd = CF * (chi - DT * (fx2 * (cen + dtaudf * dtaudf * d2taudsds / dtauds / dtauds / dtauds - d2taudfdf / dtauds) - 2. * fxsx * alfaf * AA
        -(sx2 * dtauds + 2. * fxsx * dtaudf + fx2 * dtaudf * dtaudf / dtauds) * BB))
        -DT * fx2 * CC; // chi is defined with the factor 2!!!!!!
        
    return csd;
  }
  
  void PMF::CSD(double *eta, int N1, double *csd, int N2)
  {
    double alfaf, AA, BB, CC, CF, fx2, sx2, fxsx;
    double R1, R2, R3, R4, R5, SQRTtm2, SQRTtm, SQRTtau, SQRT2t, SQRT2tau, ALFA2;
    double tm, tm2, phii, P1, cen;    
    double *EXV, *EXV2, *phi, *wa;

    assert(N1 == N2);
    EXV = (double *) malloc(N1 * sizeof(double));
    EXV2= (double *) malloc(N1 * sizeof(double));
    phi = (double *) malloc(N1 * sizeof(double));
    wa  = (double *) malloc(N1 * sizeof(double));
    
    cen = central ? 2. : 0.;
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    SQRTtau = sqrt(tau);
    SQRT2t = SQRTtau * M_SQRT2;
    SQRTtm = sqrt(tm);
    SQRTtm2 = sqrt(tm2);
    R1 = SQRTtau / SQRTtm;
    R2 = SQRT2t / SQRTtm2;
    R3 = SQRTtm / SQRTtm2;
    SQRT2tau = SQRT2t * SQRTtm2;
    ALFA2 = alfa * alfa;
    alfaf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    fx2  = cblas_ddot(3, grad_f, 1, grad_f, 1);
    sx2  = cblas_ddot(3, grad_s, 1, grad_s, 1);
    fxsx = cblas_ddot(3, grad_f, 1, grad_s, 1);
    for (int i = N1; i--;)
      EXV[i] = erfinv(eta[i]);
    for (int i = N1; i--;)
      EXV2[i] = EXV[i] * EXV[i];
    for(int i = N1; i--;)
      phi[i] = alfa + SQRT2t * EXV[i];
    for(int i = N1; i--;)
      wa[i] = exp(-EXV2[i]);
    R4 = exp(ALFA2 / 2.);
    R5 = exp(ALFA2 / 2. / tm);
    for(int i=N1; i--;)
    {
        phii = phi[i];
        BB = ALFA2 / 2. / tm /tm - (phii * EXV[i] / M_SQRT2 / SQRTtau + phii * phii / tm2 - 0.5 / tm) / tm2;
        CF = SQRTtm / SQRTtau * wa[i] * R5;
        AA = alfa / tm - phii / tm2;
        CC = M_SQRT2 / M_SQRTPI * (alfaf + phii * dtaudf / tau / tm2) * wa[i] * R4;
        csd[i] = CF * (chi - DT * (fx2 * (cen + dtaudf * dtaudf * d2taudsds / dtauds / dtauds / dtauds - d2taudfdf / dtauds) - 2. * fxsx * alfaf * AA
            -(sx2 * dtauds + 2. * fxsx * dtaudf + fx2 * dtaudf * dtaudf / dtauds) * BB))
            -DT * fx2 * CC; // chi is defined with the factor 2!!!!!!
    }
    free(EXV);
    free(EXV2);
    free(phi);
    free(wa);
  }

  double PMF::Laplace(double eta)
  {
    double R1, ALFA2, SQRT2t, phi, EXV, tm, tm2;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    ALFA2 = alfa * alfa / 2. / tm;
    R1 = - chi * M_SQRTPI * sqrt(tm) * exp(ALFA2);    
    SQRT2t = sqrt(2. * tau);
    phi = alfa + SQRT2t * EXV;    
    return R1 * (phi / tm2 + EXV / SQRT2t) * exp(- EXV * EXV / 2.);
  }

  void PMF::Laplace(double *eta, int N1, double *p, int N2)
  {
    double R1, ALFA2, SQRT2t, tm, tm2;
    double *phi,*EXV;

    assert(N1 == N2);
    EXV = (double *) malloc(N1 * sizeof(double));
    phi = (double *) malloc(N1 * sizeof(double));
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    ALFA2 = alfa * alfa / 2. / tm;
    R1 = - chi * M_SQRTPI * sqrt(tm) * exp(ALFA2);    
    SQRT2t = sqrt(2. * tau);
    for (int i = N1; i--;)
      EXV[i] = erfinv(eta[i]);
    for (int i = N1; i--;)
      phi[i] = alfa + SQRT2t * EXV[i];
    for (int i = N1; i--;)
      p[i] = R1 * (phi[i] / tm2 + EXV[i] / SQRT2t) * exp(- EXV[i] * EXV[i] / 2.);
    free(EXV);
    free(phi);
  }
  
  void PMF::CV(double eta, double *cv, int N1)
  {
    double R1, ALFA2, SQRT2t;
    double tm2, t2, alfaf;
    double phi, EXV, WA2, WA3;
    
    assert(eta >= 0 && eta <= 1);
    assert(N1 == 3);
    ALFA2 = alfa * alfa;
    EXV = erfinv(eta);
    t2 = tau + tau;
    tm2 = 1. - t2;
    SQRT2t = sqrt(t2);
    R1 = DT / tm2;
    alfaf = - M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
// Set up work arrays before computing CV. This is necessary for optimization level -O3 to work    
    phi = alfa + SQRT2t * EXV;
    WA2 = (1. + alfa * phi - phi * phi / tm2) / t2;
    WA3 = alfaf * phi;
    cv[0] = R1 * (grad_f[0] * WA3 - WA2 * (dtaudf * grad_f[0] + dtauds * grad_s[0]));
    cv[1] = R1 * (grad_f[1] * WA3 - WA2 * (dtaudf * grad_f[1] + dtauds * grad_s[1]));
    cv[2] = R1 * (grad_f[2] * WA3 - WA2 * (dtaudf * grad_f[2] + dtauds * grad_s[2]));
  }
  
  void PMF::CV(double *eta, int N1, double *cv, int N2, int N3)
  {
    double R1, ALFA2, SQRT2t;
    double tm2, t2, alfaf;
    double *phi, *EXV, *WA2, *WA3;
    
    assert(N2 == 3);
    assert(N1 == N3);
    EXV = (double *) malloc(N1 * sizeof(double));
    phi = (double *) malloc(N1 * sizeof(double));
    WA2 = (double *) malloc(N1 * sizeof(double));
    WA3 = (double *) malloc(N1 * sizeof(double));
    ALFA2 = alfa * alfa;
    for (int i = N1; i--;)
      EXV[i] = erfinv(eta[i]);
    t2 = tau + tau;
    tm2 = 1. - t2;
    SQRT2t = sqrt(t2);
    R1 = DT / tm2;
    alfaf = - M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    for (int i = N1; i--;)
      phi[i] = alfa + SQRT2t * EXV[i];
    for (int i = N1; i--;){
      WA2[i] = (1. + alfa * phi[i] - phi[i] * phi[i] / tm2) / t2;
      WA3[i] = alfaf * phi[i];}
    for (int i = N1; i--;){
      cv[i]      = R1 * (grad_f[0] * WA3[i] - WA2[i] * (dtaudf * grad_f[0] + dtauds * grad_s[0]));
      cv[i+N3]   = R1 * (grad_f[1] * WA3[i] - WA2[i] * (dtaudf * grad_f[1] + dtauds * grad_s[1]));
      cv[i+2*N3] = R1 * (grad_f[2] * WA3[i] - WA2[i] * (dtaudf * grad_f[2] + dtauds * grad_s[2]));
    }
    free(EXV);
    free(phi);
    free(WA2);
    free(WA3);
  }
}