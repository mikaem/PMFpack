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

#include "GSL_F_Root.h"

using namespace std;

namespace pmfpack{
    
  gsl_root_fsolver *
  gsl_root_fsolver_realloc (gsl_root_fsolver * s, const gsl_root_fsolver_type * T)
  {
    s->state = realloc (s->state, T->size);
    if (s->state == 0)
      {
        free (s);         /* exception in constructor, avoid memory leak */
        GSL_ERROR_VAL ("failed to reallocate space for root solver state",
                            GSL_ENOMEM, 0);
      };
    s->type = T ;
    s->function = NULL ;
    return s;
  }

  GSL_F_Root::GSL_F_Root()
  : GSLRoot()
  {
    T = gsl_root_fsolver_brent;
        
    s = gsl_root_fsolver_alloc(T);

    F.function = &func;
        
    F.params = integrator;
    
    lower = 1.e-3;
    
    upper = 1.-1.e-6;
  }

  GSL_F_Root::GSL_F_Root(double _data[], Integrator *_integrator)
  : GSLRoot(_data, _integrator)
  {
    T = gsl_root_fsolver_brent;
        
    s = gsl_root_fsolver_alloc(T);

    F.function = &func;
    
    F.params = integrator;
    
    lower = 1.e-3;
    
    upper = 1.-1.e-6;
  }

  GSL_F_Root::~GSL_F_Root()
  {
    gsl_root_fsolver_free(s);
  }

  void GSL_F_Root::realloc(const unsigned int ff)
  {
    if(ff == 1){
      T = gsl_root_fsolver_bisection;
      cout << " Using bisection bracketing solver \n" << endl;}
    else if(ff == 2){
      T = gsl_root_fsolver_brent;
      cout << " Using brent bracketing solver \n" << endl;}
    else if(ff == 3){
      T = gsl_root_fsolver_falsepos;
      cout << " Using falsepos bracketing solver \n" << endl;}
    s = gsl_root_fsolver_realloc (s, T);
    return;
  }

  double GSL_F_Root::compute(int verbose)
  {
    double x, xlo, xhi, x0;
    int status, status2, status3, iter, reset;

    (*alfa) = erfinv(1 - data[0]);
    status = gsl_root_fsolver_set (s, &F, lower, upper);
    reset = 0;
    while(status != GSL_SUCCESS && reset < 5){
      lower = lower / 10.;
      upper = 1 - (1 - upper) / 10.;
      status = gsl_root_fsolver_set (s, &F, lower, upper);
      reset++;
    }
    x = lower;
    iter = 0;
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      x0 = x;
      x = gsl_root_fsolver_root (s);
      xlo = gsl_root_fsolver_x_lower (s);
      xhi = gsl_root_fsolver_x_upper (s);
      if(xlo < 0.5){
        status2 = gsl_root_test_interval (xlo, xhi, 1e-16, 1e-15);
        status3 = gsl_root_test_delta (x, x0, 1e-16, 1e-15);
        status  = gsl_root_test_residual(data[5], 1e-15);
      }
      else{
        status2 = gsl_root_test_interval (1-xhi, 1-xlo, 1e-16, 1e-15);
        status3 = gsl_root_test_delta (1-x, 1-x0, 1e-16, 1e-15);
        status  = gsl_root_test_residual(data[5], 1e-15);
      }
      if(status2 == 0 || status3 == 0) status = 0;
              
      if (verbose){ 
        cout << "  " <<iter << "   " << x << " " << (1 - x * x) / 2 << " " << xlo << " " << xhi << endl ;
        if (status == GSL_SUCCESS)
          cout << "Converged:-) x = " << (1 - x * x) / 2 << endl;
      }
    }
    while (status == GSL_CONTINUE && iter < max_iteration);
    
    (*tau) = (1 - x * x) / 2;
    
    return (*tau);
  }
}