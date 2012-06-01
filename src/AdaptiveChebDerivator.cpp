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

#include "AdaptiveChebDerivator.h"

namespace pmfpack
{
  int pmf_cheb_init_adaptiv(gsl_cheb_series * cs, const gsl_function *func,
                    const double a, const double b, int ord, int max_ord)
  { // Should use fft here
    double sum;
    int stride;
    
    stride = max_ord / ord;

    if(a >= b) {
      GSL_ERROR_VAL("null function interval [a,b]", GSL_EDOM, 0);
    }
    cs->a = a;
    cs->b = b;
  //  cs->err = 0.0;
  
      double bma = 0.5 * (cs->b - cs->a);
      double bpa = 0.5 * (cs->b + cs->a);
      double fac = 2.0 / ord;
      
      if(cs->order_sp == cs->order)
      {
  // Set the extreme points. The central f[max_ord/2] must already have been set
          double y = 1.;
          cs->f[0] = GSL_FN_EVAL(func, (y*bma + bpa));
          y = -1.;
          cs->f[max_ord] = GSL_FN_EVAL(func, (y*bma + bpa));
          cs->order_sp = 0;
      }
      else{
        for(int k = stride; k < cs->order; k = k+2*stride) {
          double y = cos(M_PI * k / cs->order);
          cs->f[k] = GSL_FN_EVAL(func, (y*bma + bpa));
        }
      }
      
      for(int j = 0; j <= cs->order; j = j+stride) {
        sum = cs->f[0] / 2. + cs->f[cs->order] / 2. * pow(-1., j / stride);
        for(int k = stride; k < cs->order; k = k+stride) {
          sum += cs->f[k] * cos((M_PI * (j / stride) * (k*1.)) / cs->order);
                }
        cs->c[j] = fac * sum;
      }
      cs->c[0] *= 0.5;
      cs->c[cs->order] *= 0.5;
    
    return GSL_SUCCESS;
  }
  
  double pmf_cheb_eval (const gsl_cheb_series * cs, const double x, int ord, int max_ord)
  {
    double d1 = 0.0;
    double d2 = 0.0;
    double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
    double y2 = 2.0 * y;
    int stride;
    
    stride = max_ord / ord;

    for (int i = cs->order; i >= 1; i=i-stride)
      {
        double temp = d1;
        d1 = y2 * d1 - d2 + cs->c[i];
        d2 = temp;
      }

    return y * d1 - d2 + cs->c[0];
  }
  
  int pmf_cheb_calc_deriv(gsl_cheb_series * deriv, const gsl_cheb_series * f, int ord, int max_ord)
  {
    const size_t n = f->order + 1;
    const double con = 2.0 / (f->b - f->a);
    int stride;
    
    stride = max_ord / ord;
    
    if(deriv->order != f->order) 
      {
        GSL_ERROR ("order of chebyshev series must be equal", GSL_ENOMEM);
      }
    
    /* set the other parameters in the chebyshev struct */

    deriv->a = f->a;
    deriv->b = f->b;

#ifdef ERR
    deriv->err = n * n * f->c[n-1];   /* error in derivative is n^2 c_n */ 
#endif

    /* FIXME:  should probably set deriv->f[] as well */
    
    deriv->c[n-1] = 0.0;
    
    if(n > 1) {
      deriv->c[n-1-stride] = 2.0 * (n - 1.0) / (1.*stride) * f->c[n-1];

      for(int i = (n - 1 - 2 * stride); i > 0; i = i-stride) 
        deriv->c[i] = deriv->c[i + 2 * stride] + 2.0 * (i / stride + 1.0) * f->c[i+stride];

      deriv->c[0] = 0.5*(deriv->c[2*stride] + 2.0 * f->c[stride]);

      for(int i = 0; i < (n-1); i = i+stride) 
        deriv->c[i] *= con;
    }

    return GSL_SUCCESS;
  }
  
  AdaptiveChebDerivator::AdaptiveChebDerivator(Roots *_roots, int _cheb_order)
  : Derivator(_roots) 
  {
    cheb_order = pow(2, (int) log2(_cheb_order));
    assert(cheb_order == _cheb_order);
    cheb_order = _cheb_order;
    c  = gsl_cheb_alloc(cheb_order);
    c1 = gsl_cheb_alloc(cheb_order);
    c2 = gsl_cheb_alloc(cheb_order);
    F.function = &dfunction_df;
    F.params = roots;
    tol = 1.e-8;
  };
  
  AdaptiveChebDerivator::~AdaptiveChebDerivator()
  {
    gsl_cheb_free(c);
    gsl_cheb_free(c1);
    gsl_cheb_free(c2);
  }
  
  void AdaptiveChebDerivator::set_derivative_order(int order)
  {
    gsl_cheb_free(c);
    gsl_cheb_free(c1);
    gsl_cheb_free(c2);
    cheb_order = order;
    c  = gsl_cheb_alloc(cheb_order);
    c1 = gsl_cheb_alloc(cheb_order);
    c2 = gsl_cheb_alloc(cheb_order);    
  }
  
  double AdaptiveChebDerivator::compute(int verbose)
  {
    double dh, error, deriv, newerror;
    int stride;
    
    dh = gsl_min((*fmean), 1 - (*fmean));  //   0 < f < 1
    dh = gsl_min(dh, gsl_min((*sigma), (*fmean) * (1 - (*fmean)) - (*sigma)));            //   0 < sigma < f * (1 - f)
    dh = gsl_min(dh / 10., 1e-4);
    dh = pow(2, (int) log2(dh));
    
    // Compute central tau first using a robust bracketing algorithm
    roots->froot->compute(0);
    // Need to set the central tau to avoid recomputing this
    c->f[cheb_order / 2] = (*tau);
    c->order_sp = c->order;
    F.function = &dfunction_df;
    
    newerror = 1e20;
    int k = 0;
    int ord = 0;
    stride = 2;
    while(stride > 1)
    {
      ord = gsl_pow_int(2, ++k);
      stride = cheb_order / ord;
      pmf_cheb_init_adaptiv(c, &F, -dh, dh, ord, cheb_order);
      // Estimate the derivative error
      newerror = fabs(ord * ord * c->c[cheb_order - stride]);
      pmf_cheb_calc_deriv(c1, c, ord, cheb_order);   // First derivative
      deriv = pmf_cheb_eval(c1, 0, ord, cheb_order); // Its value, used for scaling the tolerance
      printf(" stride %d newerror %2.4e val %2.4e tol %2.4e\n", stride, newerror, deriv,  tol * fabs(deriv));      if (newerror < tol * fabs(deriv) || ord == cheb_order) break;
    }
    pmf_cheb_calc_deriv(c2, c1, ord, cheb_order); // Second derivative
    (*dtaudf)    = deriv;
    (*d2taudfdf) = pmf_cheb_eval(c2, 0, ord, cheb_order);
    
    c->order_sp = c->order;
    F.function = &dfunction_ds;
    newerror = 1e20;
    k = 0;
    ord = 0;
    stride = 2;
    while(stride > 1)
    {
      ord = gsl_pow_int(2, ++k);
      stride = cheb_order / ord;
      pmf_cheb_init_adaptiv(c, &F, -dh, dh, ord, cheb_order);
      // Estimate the derivative error and numerical error
      newerror = fabs(ord * ord * c->c[cheb_order - stride]);
      pmf_cheb_calc_deriv(c1, c, ord, cheb_order);   // First derivative
      deriv = pmf_cheb_eval(c1, 0, ord, cheb_order); // Its value, used for scaling the tolerance
      printf(" stride %d newerror %2.4e val %2.4e tol %2.4e\n", stride, newerror, deriv,  tol * fabs(deriv));
      if(newerror < tol * fabs(deriv) || ord == cheb_order) break;
    }
    pmf_cheb_calc_deriv(c2, c1, ord, cheb_order); // Second derivative
    (*dtauds)    = deriv;
    (*d2taudsds) = pmf_cheb_eval(c2, 0, ord, cheb_order);
  }
}
