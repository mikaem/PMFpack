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

namespace pmfpack
{  
  simple_function erfinv = &gsl_cdf_ugaussian_Pinv;
//   simple_function erfinv = &normsinv;
//   simple_function erfinv = &Erfinv;
//   simple_function erf = &Erf;  
  simple_function erf = &gsl_cdf_ugaussian_P;

  void PMF::compute_tau(int solver)
  {
    if (solver == 0)
      froot->compute(verbose);
    else if (solver == 1)
      fdfroot->compute(verbose);
    else
      std::cout << " Not implemented PMF::compute " << std::endl;
  }
  
  void PMF::reallocate_solver(Root *root, int solver)
  {
    root->realloc(solver);
  }
  
  PMF::PMF(double _fmean, double _sigma, int fsolver, int fdfsolver, int integral)
  {
    froot = NULL;
    
    fdfroot = NULL;
    
    integrator = NULL;
    
    central = true;

    verbose = 0;
    
    chi = 1.;

    set_parameters(_fmean, _sigma);
    
    set_integrator(integral);        

    set_fsolver(fsolver);
    
    set_fdfsolver(fdfsolver);
        
    fmean = &data[0];
    sigma = &data[1];
    alfa  = &data[2];    
    tau   = &data[3];
  }
  
  PMF::~PMF()
  {
    delete froot;
    delete fdfroot;
    delete integrator;
  }
  
  void PMF::set_fsolver(int fsolver)
  {
    if (froot != NULL) delete froot;    
    
    if (fsolver == 0)
    {
      froot = new GSL_F_Root(data, integrator);
    }
    else
    {
      std::cout << " Not implemented fsolver \n" << std::endl;
    }
  }
  
  void PMF::set_fdfsolver(int fdfsolver)
  {
    if (fdfroot != NULL) delete fdfroot;    
    
    if (fdfsolver == 0)
    {
      fdfroot = new GSL_FDF_Root(data, integrator);
    }
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
      integrator = new GSLIntegrator(data);
    }
//       else if (backend_integral == 1)
//       {
//         integrator = new LGIntegrator(data);
//       }

  }

  void PMF::set_parameters(double f, double s)
  {
    data[0] = f;
    data[1] = s;
    data[2] = erfinv(1 - f);
    data[3] = 0.1;
    data[4] = s + f * f;
    if(f <= 0. || f >= 1.)
    {
      std::cout << "Wrong input. 0 < fmean < 1." << std::endl;
    }
    else if(data[4] <= 0. || data[4] >= 1.){
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
  void PMF::compute_tau_and_derivatives()
  {
    double fmean0, s, alfa0, tau0, dh, dh0;
    double *sigma0;
    double tau_fp,tau_fm,tau_sp,tau_sm, s1,i1;
    int    dhi,count;

    fmean0 = data[0];
    s1 = data[1];
    i1 = data[4];
    // There is a slight difference between using central or integer moments,
    // because only the integer is involved in the computation (see intX).
    // Using central moments we change fmean and keep the central moment constant.
    // Hence the integer that is used in computations must be modified.
    // Using integer moments there's no need to do anything special, since the
    // integer moments should be kept constant when altering fmean.

    // Set the sigma0 pointer to either data[1] or data[4]
    sigma0 = central ? &data[1] : &data[4]; 
    s = (*sigma0); // remember the original value

    // Get the largest possible dh. This is computed from 0<I<1, by taking all combinations of f+-dh,s+-dh.
    // This is valid also for d2tau/dfmean/dsigma
    if (!central){
      dh = sqrt(s) - fmean0;
      dh = gsl_min(dh, 0.5 * (fmean0 - s));
      dh = gsl_min(dh, s - fmean0 * fmean0);
      dh = gsl_min(dh, ((-(2. * fmean0 - 1.) + sqrt(4. * (s - fmean0) + 1.)) / 2.));
      dh = gsl_min(dh, ((2. * fmean0 + 1. + sqrt(4. * (s + fmean0) + 1.)) / 2.));
      dh = gsl_min(dh, ((-(2. * fmean0 + 1.) + sqrt(4. * (s + fmean0) + 1.)) / 2.));   
      dh = gsl_min(dh, ((2. * fmean0 - 1. + sqrt(4. * (s - fmean0) + 1.)) / 2.));
    }
    else{
      dh = (-(2. * fmean0 - 1.) + sqrt(1. - 4. * s)) / 2.;
      dh = gsl_min(dh, (2. * fmean0 - 1. + sqrt(1. - 4. * s)) / 2);
      dh = gsl_min(dh, fmean0 - (fmean0 * fmean0 + s));
      dh = gsl_min(dh, -fmean0 + sqrt(fmean0 - s));
      dh = gsl_min(dh, fmean0 - 1. + sqrt(1. - fmean0 - s));
      dh = gsl_min(dh, -(fmean0 - 1.) + sqrt(1. - fmean0 - s));
      dh = gsl_min(dh, s);
    }
    dh = gsl_min(5.e-4, dh / 250.);
    dhi = (int) log10(dh);
    dh = pow(10, dhi);
    dh0 = dh;
        
    // Compute tau using bracketing algorithm
    compute_tau(0);
    alfa0 = data[2];
    tau0 = data[3];
    tau_fm = 1.;
    count = 0;
        
    while((tau_fm > 0.005 || tau_fm < 0.0005) && count < 10){
      if (tau_fm > 0.005) dh = dh / 10;
      else dh = dh * 10;
      data[0] = fmean0 + dh;
      dh = data[0] - fmean0; // Numerical trick
      data[1] = central ? (*sigma0) : i1 - data[0] * data[0];
      data[4] = central ? data[0] * data[0] + (*sigma0) : (*sigma0); // data[4] is modified if central moments are used
      compute_tau(1);
      tau_fp = (*tau);
      tau_fm = (fabs((*tau) - tau0) / tau0);
      count++;
    //printf("lp %f tau %f grad %f hh %2.4e\n",tau_fp,tau,tau_fm,dh);
    }
    data[3] = tau0;
    data[0] = fmean0 - dh;
    data[1] = central ? (*sigma0) : i1 - data[0] * data[0];
    data[4] = central ? data[0] * data[0] + (*sigma0) : (*sigma0);
    compute_tau(1);
    tau_fm = data[3];
    data[3] = tau0;
    data[7] = (tau_fp - tau_fm) / 2. / dh; // dtau/dfmean
    data[9] = (tau_fp - 2. * tau0 + tau_fm) / (dh * dh); //d2tau/dfmean**2     
     
    dh = dh0;
    data[0] = fmean0;
    tau_sp = 2. * tau0;
    tau_sm = (fabs(tau_sp - tau0) / tau0);
    count = 0;
    while ((tau_sm > 0.005 || tau_sm < 0.0005) && count < 10){
      if(tau_sm > 0.005) dh = dh / 10;
      else dh = dh * 10;
      (*sigma0) = s + dh;
      data[1] = central ? (*sigma0) : i1 - data[0] * data[0];
      data[4] = central ? data[0] * data[0]+(*sigma0) : (*sigma0);
      compute_tau(1);
      tau_sp = data[3];
      tau_sm = (fabs(tau_sp - tau0) / tau0);
      count++;
    }
    data[3] = tau0;
    (*sigma0) = s - dh;
    data[1] = central ? (*sigma0) : i1 - data[0] * data[0];
    data[4] = central ? data[0] * data[0] + (*sigma0) : (*sigma0);     
    compute_tau(1);
    tau_sm = data[3];
    // Reset parameters
    data[0] = fmean0;
    data[1] = s1;
    data[2] = alfa0;
    data[3] = tau0;
    data[4] = i1;
        
    // Compute derivatives
    data[8]=(tau_sp - tau_sm) / 2. / dh; // dtau/dsigma
    data[10]=(tau_sp - 2. * tau0 + tau_sm) / (dh * dh); //d2tau/dsigma**2
    
    return;
  }
  
  double PMF::PDF(double eta)
  {
    double SQRT2tau, phi, EXV;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    SQRT2tau = sqrt(2. * (*tau));
    phi = (*alfa) + SQRT2tau * EXV;
    return  SQRT2tau / sqrt(1. - 2. * (*tau)) * exp(EXV * EXV / 2. - phi * phi / (2. * (1. - 2. * (*tau))));
  }

  void PMF::PDF(double *eta, int N1, double *p, int N2)
  {
    double SQRT2tau;
    double *phi, *EXV;

    assert(N1 == N2);
    EXV = (double *) malloc(N1*sizeof(double));
    phi = (double *) malloc(N1*sizeof(double));
    SQRT2tau = sqrt(2. * (*tau));
    for (int i = 0; i < N1 ; i++){
      EXV[i] = erfinv(eta[i]);
      phi[i] = (*alfa) + SQRT2tau * EXV[i];
      p[i] = SQRT2tau / sqrt(1. - 2. * (*tau)) * exp(EXV[i] * EXV[i] / 2. - phi[i] * phi[i] / (2. * (1. - 2. * (*tau))));
    }
    free(EXV);
    free(phi);
  }
  
  double PMF::counterflow(double eta)
  {
    double tm, EXV, EXV2;

    assert(eta >= 0 && eta <= 1);
    EXV = erfinv(eta);
    EXV2 = EXV * EXV;
    tm = 1. - (*tau);
    return chi * sqrt(tm / (*tau)) * exp(-EXV2 + (*alfa) * (*alfa) / 2. / tm);
  }
  
  void PMF::counterflow(double *eta, int N1, double *p, int N2)
  {
    double tm;
    double *EXV, *phi;

    assert(N1 == N2);
    EXV = (double *) malloc(N1*sizeof(double));
    tm = 1. - (*tau);
    for (int i = 0; i < N1 ; i++){
      EXV[i] = erfinv(eta[i]);
    }
    for (int i = 0; i < N1 ; i++){
      p[i] = chi * sqrt(tm / (*tau)) * exp(-EXV[i] * EXV[i] + (*alfa) * (*alfa) / 2. / tm);
    }
    free(EXV);
  }
}