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
  
  Root::Root(double **_parameters, Integrator *_integrator)
  : integrator(_integrator), parameters(_parameters)
  {
    max_iteration = 100;
    fmean = parameters[0];
    sigma = parameters[1];
    alfa  = parameters[2];
    tau   = parameters[3];
    im    = parameters[4];
    f     = parameters[5];
    df    = parameters[6];
  }
        
  double func(double x, void *params)
  {
    double etam, etamm, etap, f1;            
    Integrator *integrator = (Integrator *) params;    
    double **data = integrator->parameters;
    
    (*data[3]) = x;
    etap = 9;
    etam = FMAX(((*data[2]) + sqrt(1 - x * x) * NIM8) / x, -8);
    etamm = ((*data[2]) + sqrt(1 - x * x) * (-NIM15)) / x;
//    etam = FMAX((data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x, -7.);
//    etam = (data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x;

    f1 = integrator->compute(&etam, &etamm, &etap);    

    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 

    return (*data[5]);
  }
  
  double dfunc(double x, void *params)
  {
    double etam, etamm, etap, f1, f2, hf;            
    Integrator *integrator = (Integrator *) params;    
    double **data = integrator->parameters;
        
    (*data[3]) = x;
    hf = FMIN(1e-5, (*data[3]) * 1e-3);
    hf = FMIN(hf, (1. - (*data[3])) * 1.e-3);
    
    etap = 9;
    etam = FMAX(((*data[2]) + sqrt(1 - x * x) * NIM8) / x, -8);
    etamm = ((*data[2]) + sqrt(1 - x * x) * (-NIM15)) / x;
//    etam = FMAX((data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x, -7.);
//    etam = (data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x;

    f1 = integrator->compute(&etam, &etamm, &etap);
    (*data[3]) -= hf;
    f2 = integrator->compute(&etam, &etamm, &etap);
    
    (*data[3]) = x;
    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 
    (*data[6]) = (f2 - f1) / hf / M_SQRTPI / M_SQRT2;

    return (*data[6]);
  }
  
  void fdfunc(double x, void *params, double *f, double *df)
  {
    double etam, etamm, etap, f1, f2, hf;            
    Integrator *integrator = (Integrator *) params;    
    double **data = integrator->parameters;
        
    (*data[3]) = x;
    hf = FMIN(1e-5, (*data[3]) * 1e-3);
    hf = FMIN(hf, (1. - (*data[3])) * 1.e-3);
    
    etap = 9;
    etam = FMAX(((*data[2]) + sqrt(1 - x * x) * NIM8) / x, -8);
    etamm = ((*data[2]) + sqrt(1 - x * x) * (-NIM15)) / x;
//    etam = FMAX((data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x, -7.);
//    etam = (data[2] + sqrt(1. - x * x) * normsinv(1.e-7)) / x;

    f1 = integrator->compute(&etam, &etamm, &etap);
    (*data[3]) -= hf;
    f2 = integrator->compute(&etam, &etamm, &etap);
    
    (*data[3]) = x;
    (*data[5]) = 1 - f1 / SQRT2PI; // Store the function value 
    (*data[6]) = (f2 - f1) / hf / SQRT2PI;
    
    (*f) = (*data[5]);
    (*df) = (*data[6]);  

    return;
  }
}