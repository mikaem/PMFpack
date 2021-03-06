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
// This is an implementation of the models described by:
//
// Mortensen and Andersson, Flow, Turbulence and Combustion Volume 76, 
//                          Number 2 (2006), 199-219)
//
// First added:  2012-05-04
// Last changed: 2012-05-04

#ifndef __PMF_H
#define __PMF_H

#include <assert.h>

#include "Boost_F_Root.h"
#include "Boost_FDF_Root.h"
#include "Root.h"
#include "GSL_F_Root.h"
#include "GSL_FDF_Root.h"
#include "Integrator.h"
#include "GSLIntegrator.h"
#include "Derivator.h"
#include "FDDerivator.h"
#include "ChebDerivator.h"
#include "AdaptiveChebDerivator.h"
#include "Lookup.h"
#include "GSLLookup.h"
#include "HDF5Lookup.h"

namespace pmfpack
{
  
  class PMF
  {
  public:
    
    PMF(double fmean, double sigma, int fsolver=0, int fdfsolver=0, int integral=0, int derivator=1);
    
    ~PMF();
    
    double fmean, sigma, alfa, tau, im, f, df;    
    double dtaudf, dtauds, d2taudfdf, d2taudsds;
    double* parameters[11];  // pointers to the parameters above
        
    double grad_f[3], grad_s[3]; // Gradients of fmean and sigma
    
    double chi;              // Mean scalar dissipation rate
    
    double DT;               // Turbulent diffusion rate
    
    Integrator *integrator;  // Class that performs the integral to determine tau (Eq. 17 of Mortensen and Andersson, FTC 2006)
    
    Roots *roots;            // struct of *froot, *fdfroot and *central
    
    Derivator *derivator;    // Class that differentiates tau wrt fmean and sigma or im
    
    int verbose;             // Print out some info if true
    
    bool central;            // True if using central moments, false if using integer moments
    
    int cheb_order;          // Order of Chebyshev approximations used in derivatives
    
    void set_parameters(double, double);
    
    void set_fmean_gradient(double, double, double dfdz=0);

    void set_sigma_gradient(double, double, double dsdz=0);
        
    void set_fsolver(int);
    
    void set_fdfsolver(int);
    
    void set_derivator(int);
    
    void reallocate_solver(Root *, int);    

    void set_integrator(int);
    
    void compute(int, bool derivatives=false);  // Compute tau with or without derivatives
    
    // Conditional models. All models are defined such that they can be computed for a single value of the mixture fraction (eta), or for an array of eta.
    double PDF(double);    
    void PDF(double *, int, double *, int);
    
    double counterflow(double);    
    void counterflow(double *, int, double *, int);
    
    double CSD(double);
    void CSD(double *, int, double *, int);
    void CSD(double *, int, double *, int, double, double, double); 
    
    double CSD_verify(double);
    void CSD_verify(double *, int, double *, int);

    double Laplace(double);
    void Laplace(double *, int, double *, int);
        
    void CV(double, double *, int);    
    void CV(double *, int, double *, int, int);
    
    double II(double);
    void II(double *, int, double *, int);
  };
  
}

#endif