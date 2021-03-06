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

#include "PMF.h"
#include "gsl/gsl_cblas.h"
#include <gsl/gsl_ieee_utils.h>

using namespace std;

namespace pmfpack
{  
  
  void PMF::reallocate_solver(Root *root, int solver)
  {
    root->realloc(solver);
  }
  
  PMF::PMF(double _fmean, double _sigma, int fsolver, int fdfsolver, int integral, int deriv)
  {
    gsl_ieee_env_setup (); /* read GSL_IEEE_MODE */    
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
    cheb_order = 8;
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
    
    // Hook up the appropriate classes to do actual calculations
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
      cout << "PMFpack:GSL fsolver" << endl; 
      roots->froot = new GSL_F_Root(integrator);
    }
#ifdef HAS_BOOST    
    else if (fsolver == 1)
    {
      cout << "PMFpack:Boost fsolver" << endl; 
      roots->froot = new Boost_F_Root(integrator);
    }
#endif
    else
    {
      cout << "PMFpack:Not implemented fsolver \n" << endl;
    }
  }
  
  void PMF::set_fdfsolver(int fdfsolver)
  {
    if (roots->fdfroot != NULL) delete roots->fdfroot;    
    
    if (fdfsolver == 0)
    {
      cout << "PMFpack:GSL fdfsolver" << endl;
      roots->fdfroot = new GSL_FDF_Root(integrator);
    }
#ifdef HAS_BOOST
    else if (fdfsolver == 1)
    {
      cout << "PMFpack:Boost fdfsolver" << endl; 
      roots->fdfroot = new Boost_FDF_Root(integrator);
    }
#endif
    else
    {
      cout << "PMFpack:Not implemented fdfsolver \n" << endl;
    }
  }
  
  void PMF::set_integrator(int integral)
  {
    if (integrator != NULL) delete integrator;
    
    if (integral == 0)
    {
      cout << "PMFpack:GSL integrator" << endl; 
      integrator = new GSLIntegrator(parameters);
    }
    else
    {
      cout << "PMFpack:Not implemented integrator \n" << endl;
    }
  }
  
  void PMF::set_derivator(int deriv)
  {
    if (derivator != NULL) delete derivator;
    
    if (deriv == 0)
    {
      cout << "PMFpack:Finite difference derivator" << endl; 
      derivator = new FDDerivator(roots);
    }
    else if (deriv == 1)
    {
      cout << "PMFpack:Chebyshev derivator" << endl; 
      derivator = new ChebDerivator(roots, cheb_order);
    }
    else if (deriv == 2)
    {
      cout << "PMFpack:Adaptive Chebyshev derivator" << endl; 
      derivator = new AdaptiveChebDerivator(roots, 16); // Maximum order 16
    }
    else
    {
      cout << "PMFpack:Not implemented derivator \n" << endl;
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
      cout << "PMFpack:Wrong input. 0 < fmean = (" << f << ") < 1" << endl;
    }
    else if(s <= 0. || s / f / (1 - f) >= 1.)
    {
      cout << "PMFpack:Wrong input. 0 < Intensity of segregation (" << s / f / (1-f) << ") < 1" << endl;
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
        cout << " Not implemented PMF::compute " << endl;
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
    double dalfadf, AA, BB, CC, CF, fx2, sx2, fxsx;
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
    dalfadf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
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
    CC = M_SQRT2 / M_SQRTPI * (dalfadf + phi * dtaudf / tau / tm2) * wa * R4;
    csd = CF * (chi - DT * (fx2 * (cen + dtaudf * dtaudf * d2taudsds / dtauds / dtauds / dtauds - d2taudfdf / dtauds) - 2. * fxsx * dalfadf * AA
        -(sx2 * dtauds + 2. * fxsx * dtaudf + fx2 * dtaudf * dtaudf / dtauds) * BB))
        -DT * fx2 * CC; // chi is defined with the factor 2!!!!!!
        
    return csd;
  }
  
  void PMF::CSD(double *eta, int N1, double *csd, int N2)
  {
    double dalfadf, AA, BB, CC, CF, fx2, sx2, fxsx;
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
    dalfadf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
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
        CC = M_SQRT2 / M_SQRTPI * (dalfadf + phii * dtaudf / tau / tm2) * wa[i] * R4;
        csd[i] = CF * (chi - DT * (fx2 * (cen + dtaudf * dtaudf * d2taudsds / dtauds / dtauds / dtauds - d2taudfdf / dtauds) - 2. * fxsx * dalfadf * AA
            -(sx2 * dtauds + 2. * fxsx * dtaudf + fx2 * dtaudf * dtaudf / dtauds) * BB))
            -DT * fx2 * CC; // chi is defined with the factor 2!!!!!!
    }
    free(EXV);
    free(EXV2);
    free(phi);
    free(wa);
  }
  
  void PMF::CSD(double *eta, int N1, double *csd, int N2, double fx2, double sx2, double fxsx)
  {
    // With this routine you can provide the dot products of the gradients directly
    double dalfadf, AA, BB, CC, CF;
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
    dalfadf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
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
        CC = M_SQRT2 / M_SQRTPI * (dalfadf + phii * dtaudf / tau / tm2) * wa[i] * R4;
        csd[i] = CF * (chi - DT * (fx2 * (cen + dtaudf * dtaudf * d2taudsds / dtauds / dtauds / dtauds - d2taudfdf / dtauds) - 2. * fxsx * dalfadf * AA
            -(sx2 * dtauds + 2. * fxsx * dtaudf + fx2 * dtaudf * dtaudf / dtauds) * BB))
            -DT * fx2 * CC; // chi is defined with the factor 2!!!!!!
    }
    free(EXV);
    free(EXV2);
    free(phi);
    free(wa);      
  }
  
  double PMF::CSD_verify(double eta)
  {
    double dalfadf, fx2, sx2, fxsx;
    double SQRTtm2, SQRTtm, SQRTtau, SQRT2t, ALFA2;
    double tm, tm2, cen, pdf, N, dpdfdtau, dNdtau, dpdfda, dNda;    
    double EXV, EXV2, phi, csd, ee;
    double dIIdtau, d2IIdtau2, dIIdalfa, d2IIdadtau, dIIds;
    double d2IIdsds, d2IIdfdf, d2IIdfds;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    cen = central ? 2. : 0.;
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    SQRTtau = sqrt(tau);
    SQRT2t = SQRTtau * M_SQRT2;
    SQRTtm = sqrt(tm);
    SQRTtm2 = sqrt(tm2);
    ALFA2 = alfa * alfa;
    dalfadf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    fx2  = cblas_ddot(3, grad_f, 1, grad_f, 1);
    sx2  = cblas_ddot(3, grad_s, 1, grad_s, 1);
    fxsx = cblas_ddot(3, grad_f, 1, grad_s, 1);
    phi = alfa + SQRT2t * EXV;
    EXV2 = EXV * EXV;
    ee = (EXV / M_SQRT2 + SQRTtau * alfa) / SQRTtm2;
    pdf = SQRT2t / SQRTtm2 * exp(EXV2 / 2. - phi * phi / (2. * tm2));
    N = SQRTtm / SQRTtau * exp(-EXV2 + ALFA2 / 2. / tm);
    dpdfda = - phi / tm2 * pdf;
    dpdfdtau = pdf / 2. / tau / tm2 * (1. + phi * alfa - phi * phi / tm2);
    dNda = alfa / tm * N;
    dNdtau = N * (ALFA2 / 2. / tm / tm - 1. / 2. / tau / tm);
    dIIds = 0.5 * pdf * N;
    dIIdtau = dIIds / dtauds;
    d2IIdtau2 = 0.5 * (pdf * dNdtau + N * dpdfdtau) / dtauds - dIIdtau * d2taudsds / pow(dtauds, 2);
    dIIdalfa = - 1. / dalfadf * erf((EXV / M_SQRT2 + SQRTtau * alfa) * M_SQRT2 / SQRTtm2);    
    d2IIdadtau = - phi / 2. / dalfadf / sqrt(tau * M_PI) / tm2 / SQRTtm2 * exp(-ee * ee);
    d2IIdsds = 0.5 * (pdf * dNdtau + N * dpdfdtau) * dtauds;
    d2IIdfds = 0.5 * (pdf * (dNda * dalfadf + dNdtau * dtaudf) 
                      + N * (dpdfda * dalfadf + dpdfdtau * dtaudf));
    d2IIdfdf = - dalfadf * sqrt(tau / M_PI / tm2) * exp(- ee * ee) 
               + d2IIdtau2 * dtaudf * dtaudf + dIIdtau * d2taudfdf 
               + 2. * d2IIdadtau * dtaudf * dalfadf;
              
    csd =  N * (chi - DT * fx2 * cen) 
        + 2. * DT / pdf * (
                d2IIdfdf * fx2
              + d2IIdsds * sx2
         + 2. * d2IIdfds * fxsx); // chi is defined with the factor 2!!!!!!
        
    return csd;
  }
  
  void PMF::CSD_verify(double *eta, int N1, double *csd, int N2)
  {
    // This is a more detailed version of CSD, just to remind me what is actually computed (MM).
    double dalfadf, fx2, sx2, fxsx;
    double SQRTtm2, SQRTtm, SQRTtau, SQRT2t, ALFA2;
    double tm, tm2, cen, dpdfda, dNda, dIIds;    
    double *dIIdtau, *d2IIdtau2, *dIIdalfa, *d2IIdadtau;
    double *d2IIdsds, *d2IIdfdf, *d2IIdfds;
    double *EXV, *EXV2, *phi, *wa, *pdf, *N, *dpdfdtau;
    double *dNdtau;
    
    assert(N1 == N2);
    EXV = (double *) malloc(15 * N1 * sizeof(double));
    EXV2= EXV + N1;
    phi = EXV + 2 * N1;
    wa  = EXV + 3 * N1;
    pdf = EXV + 4 * N1;
    N   = EXV + 5 * N1;
    dpdfdtau   = EXV + 6 * N1;
    dNdtau     = EXV + 7 * N1;
    dIIdtau    = EXV + 8 * N1;
    d2IIdtau2  = EXV + 9 * N1;
    dIIdalfa   = EXV + 10 * N1;
    d2IIdadtau = EXV + 11 * N1;
    d2IIdsds   = EXV + 12 * N1;
    d2IIdfdf   = EXV + 13 * N1;
    d2IIdfds   = EXV + 14 * N1;
            
    cen = central ? 2. : 0.;
    tm = 1. - tau;
    tm2 = 1. - 2. * tau;
    SQRTtau = sqrt(tau);
    SQRT2t = SQRTtau * M_SQRT2;
    SQRTtm = sqrt(tm);
    SQRTtm2 = sqrt(tm2);
    ALFA2 = alfa * alfa;
    dalfadf = -M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    fx2  = cblas_ddot(3, grad_f, 1, grad_f, 1);
    sx2  = cblas_ddot(3, grad_s, 1, grad_s, 1);
    fxsx = cblas_ddot(3, grad_f, 1, grad_s, 1);
    
    for (int i = N1; i--;)
    {
      EXV[i] = erfinv(eta[i]);
      EXV2[i] = EXV[i] * EXV[i];
      phi[i] = alfa + SQRT2t * EXV[i];
      wa[i] = (EXV[i] / M_SQRT2 + SQRTtau * alfa) / SQRTtm2;
      pdf[i] = SQRT2t / SQRTtm2 * exp(EXV2[i] / 2. - phi[i] * phi[i] / (2. * tm2));
      N[i] = SQRTtm / SQRTtau * exp(-EXV2[i] + ALFA2 / 2. / tm);
      dpdfdtau[i] = pdf[i] / 2. / tau / tm2 * (1. + phi[i] * alfa - phi[i] * phi[i] / tm2);
      dNdtau[i] = N[i] * (ALFA2 / 2. / tm / tm - 1. / 2. / tau / tm);
      dIIds = 0.5 * pdf[i] * N[i];
      dIIdtau[i] = dIIds / dtauds;
    }
    for(int i = N1; i--;)
    {
      d2IIdtau2[i] = 0.5 * (pdf[i] * dNdtau[i] + N[i] * dpdfdtau[i]) / dtauds - dIIdtau[i] * d2taudsds / pow(dtauds, 2);
      dIIdalfa[i] = - 1. / dalfadf * erf(wa[i] * M_SQRT2);
      d2IIdadtau[i] = - phi[i] / 2. / dalfadf / sqrt(tau * M_PI) / tm2 / SQRTtm2 * exp(-wa[i] * wa[i]);
    }
    for(int i = N1; i--;)
    {
      dNda = alfa / tm * N[i];
      dpdfda = - phi[i] / tm2 * pdf[i];
      d2IIdsds[i] = 0.5 * (pdf[i] * dNdtau[i] + N[i] * dpdfdtau[i]) * dtauds;
      d2IIdfds[i] = 0.5 * (pdf[i] * (dNda * dalfadf + dNdtau[i] * dtaudf) 
                      + N[i] * (dpdfda * dalfadf + dpdfdtau[i] * dtaudf));
      d2IIdfdf[i] = - dalfadf * sqrt(tau / M_PI / tm2) * exp(- wa[i] * wa[i]) 
               + d2IIdtau2[i] * dtaudf * dtaudf + dIIdtau[i] * d2taudfdf 
               + 2. * d2IIdadtau[i] * dtaudf * dalfadf;
    }
    for(int i = N1; i--;)
    {
      csd[i] =  N[i] * (chi - DT * fx2 * cen) 
          + 2. * DT / pdf[i] * (
                  d2IIdfdf[i] * fx2
                + d2IIdsds[i] * sx2
           + 2. * d2IIdfds[i] * fxsx); // chi is defined with the factor 2!!!!!!
    }
    
    free(EXV);
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
    double tm2, t2, dalfadf;
    double phi, EXV, WA2, WA3;
    
    assert(eta >= 0 && eta <= 1);
    assert(N1 == 3);
    ALFA2 = alfa * alfa;
    EXV = erfinv(eta);
    t2 = tau + tau;
    tm2 = 1. - t2;
    SQRT2t = sqrt(t2);
    R1 = DT / tm2;
    dalfadf = - M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    phi = alfa + SQRT2t * EXV;
    WA2 = (1. + alfa * phi - phi * phi / tm2) / t2;
    WA3 = dalfadf * phi;
    cv[0] = R1 * (grad_f[0] * WA3 - WA2 * (dtaudf * grad_f[0] + dtauds * grad_s[0]));
    cv[1] = R1 * (grad_f[1] * WA3 - WA2 * (dtaudf * grad_f[1] + dtauds * grad_s[1]));
    cv[2] = R1 * (grad_f[2] * WA3 - WA2 * (dtaudf * grad_f[2] + dtauds * grad_s[2]));
  }
  
  void PMF::CV(double *eta, int N1, double *cv, int N2, int N3)
  {
    double R1, ALFA2, SQRT2t;
    double tm2, t2, dalfadf;
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
    dalfadf = - M_SQRT2 * M_SQRTPI * exp(ALFA2 / 2.);
    for (int i = N1; i--;)
      phi[i] = alfa + SQRT2t * EXV[i];
    for (int i = N1; i--;){
      WA2[i] = (1. + alfa * phi[i] - phi[i] * phi[i] / tm2) / t2;
      WA3[i] = dalfadf * phi[i];}
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
  
  double xr(double phi, void * par)
  {
    double **data = (double **) par;
    double tau = (*data[3]);
    double alfa = (*data[2]);
    return erf((phi - alfa) / sqrt(2. * tau)) / sqrt(2. * M_PI * (1. - 2. * tau)) * exp(- phi * phi / 2. / (1. - 2. * tau));
  }
  
  double PMF::II(double eta)
  {
    // Return the double integral II(eta) = \int_{0}^{eta} \int_{0}^{eta'} pdf(eta'') deta'' deta'
    // Integration by parts give II(eta) = eta * \int_{0}^{eta} pdf(eta') deta' - \int_{0}^{eta} X * pdf * deta'
    // Transformation to mapped space: II(eta) = eta * \int_{-\infty}^{phi} r(phi') dphi' - \int_{-\infty}^{phi} X(phi') * r(phi') dphi'
    // Or simply as:
    //
    //   II(eta) = eta * erf(phi / sqrt(1. - 2. * tau)) - \int_{-\infty}^{phi} X(phi') * r(phi') dphi' 
    //
    // where erf as before really is the normal distribution (see PMFMath.cpp) 
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(96);
    gsl_function F;
    
    F.function = &xr;
    F.params = parameters;
    double phi = alfa + sqrt(2. * tau) * erfinv(eta);
    double integral = eta * erf(phi / sqrt(1. - 2. * tau));
    double phi_min = alfa - 12. * sqrt(2. * tau);
    integral -= gsl_integration_glfixed (&F, phi_min, phi, table);
    gsl_integration_glfixed_table_free(table);
    return integral;
  }
  
  void PMF::II(double *eta, int N1, double *x, int N2)
  {
    // Return the double integral II(eta) = \int_{0}^{eta} \int_{0}^{eta'} pdf(eta'') deta'' deta'
    // Integration by parts give II(eta) = eta * \int_{0}^{eta} pdf(eta') deta' - \int_{0}^{eta} X * pdf * deta'
    // Transformation to mapped space: II(eta) = eta * \int_{-\infty}^{phi} r(phi') dphi' - \int_{-\infty}^{phi} X(phi') * r(phi') dphi'
    // Or simply as:
    //
    //   II(eta) = eta * erf(phi / sqrt(1. - 2. * tau)) - \int_{-\infty}^{phi} X(phi') * r(phi') dphi' 
    //
    // where erf as before really is the normal distribution (see PMFMath.cpp) 
    double phi, phi_min;
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(96);
    gsl_function F;
    
    F.function = &xr;
    F.params = parameters;
    
    phi_min = alfa - 12. * sqrt(2. * tau);
    for (int i=0; i<N1; i++)
    {
      phi = alfa + sqrt(2. * tau) * erfinv(eta[i]);
      x[i] = eta[i] * erf(phi / sqrt(1. - 2. * tau));
      x[i] -= gsl_integration_glfixed (&F, phi_min, phi, table);
    }    
    gsl_integration_glfixed_table_free(table);
  }
}
