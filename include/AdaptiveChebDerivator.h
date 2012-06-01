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
// First added:  2012-06-01
// Last changed: 2012-06-01
//
// This class reimplements GSL's Chebyshev routines to help speed
// things up and avoid recomputing the central tau.
// We allocate the gsl_cheb_serieses for a maximum of cheb_order,
// that must be an integer power of 2^n (4, 8, 16, 32, ..)
// We then start by using 4'th order and compute the error
// If the error is too large we move to 8'th order - reusing 
// all previously computed function values. Compute and check 
// the error and so on.

#ifndef __ADAPTIVECHEBDERIVATOR_H
#define __ADAPTIVECHEBDERIVATOR_H

#include "Derivator.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>

namespace pmfpack
{
  int    pmf_cheb_calc_deriv(gsl_cheb_series *, const gsl_cheb_series *, int, int);
  double pmf_cheb_eval (const gsl_cheb_series *, const double, int, int);
  int    pmf_cheb_init_adaptiv(gsl_cheb_series *, const gsl_function *,
                    const double, const double, int, int);
  
  class AdaptiveChebDerivator : public Derivator 
  {
  public:
    
    AdaptiveChebDerivator(Roots *, int);
            
    ~AdaptiveChebDerivator();
      
    virtual double compute(int);
    
    void set_derivative_order(int);
    
    int cheb_order;
    
    double tol;
    
    gsl_cheb_series *c, *c1, *c2;
    
    gsl_function F;

  };
}

#endif
