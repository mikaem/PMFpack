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

#include <iostream>
#include <assert.h>

#include "Root.h"
#include "GSL_F_Root.h"
#include "GSL_FDF_Root.h"
#include "Integrator.h"
#include "GSLIntegrator.h"

namespace pmfpack
{

  class PMF
  {
  public:
    PMF() {};
    
    PMF(double fmean, double sigma, int fsolver=0, int fdfsolver=0, int integral=0);
    
    ~PMF();
    
    // [fmean, sigma, alfa, tau, Is, f, df, dtaudf, dtauds, d2taudf2, d2tauds2, d2taudfds]
    double data[12];
    
    double *fmean, *sigma, *alfa, *tau;
    
    double grad_f[3], grad_s[3]; // Gradients of fmean and sigma
    
    double chi;    // Mean scalar dissipation rate
    
    double DT;     // Turbulent diffusion rate
    
    Integrator *integrator;  // Class that performs the integral to determine tau (Eq. 17 of Mortensen and Andersson, FTC 2006)
    
    Root *froot, *fdfroot;
    
    int verbose;
    
    bool central;
    
    void set_parameters(double, double);
    
    void set_fmean_gradient(double, double, double dfdz=0);

    void set_sigma_gradient(double, double, double dsdz=0);
        
    void set_fsolver(int);
    
    void set_fdfsolver(int);
    
    void reallocate_solver(Root *, int);    

    void set_integrator(int);
    
    void compute_tau(int);
    
    void compute_tau_and_derivatives();
          
    double PDF(double);
    
    void PDF(double *, int, double *, int);
    
    double counterflow(double);
    
    void counterflow(double *, int, double *, int);
    
  };
}

#endif