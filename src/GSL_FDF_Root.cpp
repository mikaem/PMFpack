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

#include "GSL_FDF_Root.h"

using namespace std;

namespace pmfpack
{

  gsl_root_fdfsolver *
  gsl_root_fdfsolver_realloc (gsl_root_fdfsolver * s, const gsl_root_fdfsolver_type * T)
  {
    s->state = realloc (s->state, T->size);
    if (s->state == 0)
      {
        free (s);         /* exception in constructor, avoid memory leak */
        GSL_ERROR_VAL ("failed to allocate space for root solver state",
                          GSL_ENOMEM, 0);
      };
    s->type = T ;
    s->fdf = NULL;
    return s;
  }

  GSL_FDF_Root::GSL_FDF_Root()
  : GSLRoot()
  {
      T = gsl_root_fdfsolver_newton;      
      s = gsl_root_fdfsolver_alloc(T);      
      F.f = &gslfunc;      
      F.df = &gsldfunc;            
      F.fdf = &gslfdfunc;      
      F.params = integrator;    
  }
  
  GSL_FDF_Root::GSL_FDF_Root(Integrator *_integrator)
  : GSLRoot(_integrator)
  {
      T = gsl_root_fdfsolver_newton;      
      s = gsl_root_fdfsolver_alloc(T);      
      F.f = &gslfunc;      
      F.df = &gsldfunc;            
      F.fdf = &gslfdfunc;      
      F.params = integrator;    
  }
  
  GSL_FDF_Root::~GSL_FDF_Root()
  {
    gsl_root_fdfsolver_free(s);
  }

  double GSL_FDF_Root::compute(int verbose)
  {
    int status, status2, status3, iter;
    double x, x0, xhi, Is;
    double tau0;
    
    status = 0;
    iter = 0;
    xhi = 1;
    tau0 = (*tau);

// Remember we are solving for sr = sqrt(1-2*tau)
// The initial guess is for tau. Hence this must first be converted to sr
    (*alfa) = erfinv(1 - (*fmean));
    (*tau)  = sqrt(1 - 2 * (*tau));
    x0 = x = (*tau); // Initial guess    
    gsl_root_fdfsolver_set (s, &F, x0);
    status2 = gsl_root_test_residual((*f), 1e-15);
    status = (fabs((*f) / (*df)) < 1e-16) ? 0 : 1;
    if (status == 0 || status2 == 0){
      if (verbose == 1){
        printf("Initial value accepted (f = %2.4e). \n", *f);
      }
      (*tau) = tau0; // Transform back to tau
      return -1;
    }
    if(verbose == 1){
      std::cout << "Root using " << gsl_root_fdfsolver_name(s) << " method" << std::endl;
      std::cout << "Iter      root         error             f  " << std::endl;
    }
    do  
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      if (x < 0.5){
        status = gsl_root_test_delta (x, x0, 1e-15, 1e-14);
        status2 = gsl_root_test_residual((*f), 1e-15);
        status3 = (fabs((*f) / (*df)) < 1e-16) ? 0 : 1;}
      else{
        status = gsl_root_test_delta (1-x, 1-x0, 1e-15, 1e-14);
        status2 = gsl_root_test_residual((*f), 1e-15);
        status3 = (fabs((*f) / (*df)) < 1e-16) ? 0 : 1;}
      if (status == GSL_SUCCESS || status2 == GSL_SUCCESS || status3 == GSL_SUCCESS)
        status = GSL_SUCCESS;
      
      if (verbose == 1){
        std::cout << setprecision(6) << std::setw(4) << iter << std::setw(14) << 0.5*(1.-x*x) << std::setw(14) << x-x0 << std::setw(14) << (*f) << std::endl;
      if (status == GSL_SUCCESS)
        std::cout << "Converged, root = " << 0.5 * (1. - x * x) << std::endl;
      }
    }
    while (status == GSL_CONTINUE && iter < max_iteration && x > 0. && x < xhi);

    if(status != GSL_SUCCESS){
      printf(" Error, status = %d %d %d\n", status, status2, status3);
      return -1;
    }
    (*tau) = (1 - x * x) / 2; // Transform back to tau
    return (*tau);
  }
  
  void GSL_FDF_Root::realloc(const unsigned int fdf)
  {
    if (fdf == 1){
      T = gsl_root_fdfsolver_newton;
      cout << " Using Newton FDF solver \n" << endl;}
    else if (fdf == 2){
      T = gsl_root_fdfsolver_secant;
      cout << " Using Secant FDF solver \n" << endl;}
    else if (fdf == 3){
      T = gsl_root_fdfsolver_steffenson;
      cout << " Using Steffenson FDF solver \n" << endl;}
    s = gsl_root_fdfsolver_realloc (s, T);
    return;
  }
}