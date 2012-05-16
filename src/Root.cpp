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

#include "Root.h"

using namespace std;

namespace pmfpack
{
    
  Root::Root()
  {
    max_iteration = 100;    
  }
  
  Root::Root(Integrator *_integrator)
  : integrator(_integrator)
  {
    parameters = integrator->parameters;
    max_iteration = 100;
    fmean = parameters[0];
    sigma = parameters[1];
    alfa  = parameters[2];
    tau   = parameters[3];
    im    = parameters[4];
    f     = parameters[5];
    df    = parameters[6];    
  }
        
  double func(double x, Integrator *integrator)
  {
    double **data = integrator->parameters;
    
    (*data[3]) = x;

    double f1 = integrator->compute();

    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 

    return (*data[5]);
  }
  
  double dfunc(double x, Integrator *integrator)
  {
    double f1, f2, hf;            
    double **data = integrator->parameters;
        
    (*data[3]) = x;
    hf = FMIN(1e-5, (*data[3]) * 1e-3);
    hf = FMIN(hf, (1. - (*data[3])) * 1.e-3);
    
    f1 = integrator->compute();
    (*data[3]) -= hf;
    f2 = integrator->compute();
    
    (*data[3]) = x;
    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 
    (*data[6]) = (f2 - f1) / hf / M_SQRTPI / M_SQRT2;

    return (*data[6]);
  }
  
  void fdfunc(double x, Integrator *integrator, double *f, double *df)
  {
    double f1, f2, hf;            
    double **data = integrator->parameters;
        
    (*data[3]) = x;
    hf = FMIN(1e-5, (*data[3]) * 1e-3);
    hf = FMIN(hf, (1. - (*data[3])) * 1.e-3);
    
    f1 = integrator->compute();
    (*data[3]) -= hf;
    f2 = integrator->compute();
    
    (*data[3]) = x;
    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 
    (*data[6]) = (f2 - f1) / hf / SQRT2PI;
    
    (*f) = (*data[5]);
    (*df) = (*data[6]);  
  }
  
  void fd2func(double x, Integrator *integrator, double *f, double *df, double *d2f)
  {
    double f1, f2, f3, hf;            
    double **data = integrator->parameters;
        
    (*data[3]) = x;
    hf = FMIN(1e-5, (*data[3]) * 1e-3);
    hf = FMIN(hf, (1. - (*data[3])) * 1.e-3);
    
    f1 = integrator->compute();
    (*data[3]) -= hf;
    f2 = integrator->compute();
    (*data[3]) += 2 * hf;
    f3 = integrator->compute();
    
    (*data[3]) = x;
    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 
    (*data[6]) = (f2 - f3) / 2. / hf / SQRT2PI;    
    (*d2f) = - (f2 - 2 * f1 + f3) / hf / hf / SQRT2PI;    
    (*f) = (*data[5]);
    (*df) = (*data[6]);  
  }
}