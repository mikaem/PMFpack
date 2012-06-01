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
// First added:  2012-05-14
// Last changed: 2012-05-14

#include "ChebDerivator.h"

namespace pmfpack
{
  ChebDerivator::ChebDerivator(Roots *_roots, int _cheb_order)
  : Derivator(_roots) 
  {
    cheb_order = _cheb_order;
    c  = gsl_cheb_alloc(cheb_order);
    c1 = gsl_cheb_alloc(cheb_order);
    c2 = gsl_cheb_alloc(cheb_order);
    F.function = &dfunction_df;
    F.params = roots;
  };
  
  ChebDerivator::~ChebDerivator()
  {
    gsl_cheb_free(c);
    gsl_cheb_free(c1);
    gsl_cheb_free(c2);
  }
  
  void ChebDerivator::set_derivative_order(int order)
  {
    gsl_cheb_free(c);
    gsl_cheb_free(c1);
    gsl_cheb_free(c2);
    cheb_order = order;
    c  = gsl_cheb_alloc(cheb_order);
    c1 = gsl_cheb_alloc(cheb_order);
    c2 = gsl_cheb_alloc(cheb_order);    
  }
  
  double ChebDerivator::compute(int verbose)
  {
    double dh, aerr, error, error1;
    
    dh = gsl_min((*fmean), 1 - (*fmean));  //   0 < f < 1
    dh = gsl_min(dh, gsl_min((*sigma), (*fmean) * (1 - (*fmean)) - (*sigma)));            //   0 < sigma < f * (1 - f)
    dh = gsl_min(dh / 10., 1e-4);
    dh = pow(2, (int) log2(dh));
    
    // Compute central tau first using a robust bracketing algorithm
    roots->froot->compute(0);
    F.function = &dfunction_df;
    gsl_cheb_init(c, &F, -dh, dh);
    gsl_cheb_calc_deriv(c1, c);  // First derivative
    gsl_cheb_calc_deriv(c2, c1); // Second derivative
    if (verbose){
      error = fabs(cheb_order * cheb_order * c->c[cheb_order-1]); // Approx error in first derivative
      //error1 = fabs(cheb_order * cheb_order * c1->c[cheb_order-1]); // Approx error in second derivative
      aerr = 0.;
      for(int i=0; i<=cheb_order; i++) aerr += fabs(c->c[i]);
      std::cout << "Chebyshev error fmean: derivative = " << error <<  ", Numerical error " << aerr*GSL_DBL_EPSILON << std::endl;
    }
    (*dtaudf)    = gsl_cheb_eval(c1, 0);
    (*d2taudfdf) = gsl_cheb_eval(c2, 0);
    
    F.function = &dfunction_ds;
    gsl_cheb_init(c, &F, -dh, dh);    
    gsl_cheb_calc_deriv(c1, c);
    gsl_cheb_calc_deriv(c2, c1);
    if (verbose){
      error = fabs(cheb_order * cheb_order * c->c[cheb_order-1]); // Approx error in derivative
      //error1 = fabs(cheb_order * cheb_order * c1->c[cheb_order-1]); // Approx error in second derivative
      aerr = 0.;
      for(int i=0; i<=cheb_order; i++) aerr += fabs(c1->c[i]);
      std::cout << "Chebyshev error sigma: derivative = " << error <<  ", Numerical error " << aerr*GSL_DBL_EPSILON << std::endl;
    }
    (*dtauds)    = gsl_cheb_eval(c1, 0);
    (*d2taudsds) = gsl_cheb_eval(c2, 0);    
  }
}
