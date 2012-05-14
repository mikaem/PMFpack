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

#include "ChebDerivator.h"
#include <iostream>

namespace pmfpack
{
  ChebDerivator::ChebDerivator(bool _central, Root *_froot, Root *_fdfroot, int _cheb_order)
  : Derivator(_central, _froot, _fdfroot) 
  {
    cheb_order = _cheb_order;
    c  = gsl_cheb_alloc(cheb_order);
    c1 = gsl_cheb_alloc(cheb_order);
    c2 = gsl_cheb_alloc(cheb_order);
    F.function = &dfunction_df;
    F.params = params;
  };
  
  ChebDerivator::~ChebDerivator()
  {
    gsl_cheb_free(c);
    gsl_cheb_free(c1);
    gsl_cheb_free(c2);
  }
  
  double ChebDerivator::compute(int verbose)
  {
    double fmean0, sigma0, alfa0, tau0, im0, dh, dh0, aerr;
    double tau_p, tau_m, s, s1,i1, error;
    int    dhi, count;

    fmean0 = (*fmean);
    im0 = (*im);
    sigma0 = (*sigma);    
//     fmean0 = gsl_min(fmean0, 1-fmean0);
    dh = gsl_min(sigma0, fmean0 - im0);
    dh = gsl_min(dh, fabs(((1 - 2 * fmean0) + sqrt(1 - 4 * sigma0)) / 2));
    dh = gsl_min(dh, fabs(((1 - 2 * fmean0) - sqrt(1 - 4 * sigma0)) / 2));
    dh = gsl_min(1e-3, dh * 0.2);
    
    // Compute central tau first using a robust bracketing algorithm
    params->froot->compute(0);
    F.function = &dfunction_df;
    gsl_cheb_init(c, &F, -dh, dh);
    gsl_cheb_calc_deriv(c1, c);  // First derivative
    gsl_cheb_calc_deriv(c2, c1); // Second derivative
    error = fabs(cheb_order * cheb_order * c->c[cheb_order-1]); // Approx error in first derivative
//    error1 = fabs(cheb_order * cheb_order * c1->c[cheb_order-1]); // Approx error in second derivative
    if (verbose){
      aerr = 0.;
      for(int i=0; i<=cheb_order; i++) aerr += fabs(c->c[i]);
      std::cout << "Chebyshev error fmean: derivative = " << error <<  ", Numerical error " << aerr*GSL_DBL_EPSILON << std::endl;
    }
    (*tau) = gsl_cheb_eval(c, 0);
    (*dtaudf)    = gsl_cheb_eval(c1, 0);
    (*d2taudfdf) = gsl_cheb_eval(c2, 0);
    
    F.function = &dfunction_ds;
    gsl_cheb_init(c, &F, -dh, dh);    
    gsl_cheb_calc_deriv(c1, c);
    gsl_cheb_calc_deriv(c2, c1);
    error = fabs(cheb_order * cheb_order * c->c[cheb_order-1]); // Approx error in derivative
//     error1 = fabs(cheb_order * cheb_order * c1->c[cheb_order-1]); // Approx error in second derivative
    if (verbose){
      aerr = 0.;
      for(int i=0; i<=cheb_order; i++) aerr += fabs(c1->c[i]);
      std::cout << "Chebyshev error sigma: derivative = " << error <<  ", Numerical error " << aerr*GSL_DBL_EPSILON << std::endl;
    }
    (*dtauds)    = gsl_cheb_eval(c1, 0);
    (*d2taudsds) = gsl_cheb_eval(c2, 0);
    
  }
}
