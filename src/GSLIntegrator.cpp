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

#include "GSLIntegrator.h"
#include <iostream>

namespace pmfpack
{

  GSLIntegrator::GSLIntegrator(double **_parameters)
  : Integrator(_parameters)
  {
    F.function = &function;    
    F.params = parameters;
    w = gsl_integration_workspace_alloc (1000);        
  }
  
  GSLIntegrator::~GSLIntegrator()
  {
    gsl_integration_workspace_free(w);
  }
  
  double GSLIntegrator::compute()
  {
    int status;
    double etam, etamm, etap, f1, f2, aerr, error;    
    size_t neval;
    
    aerr = 5.e-8 / (SQRT2PI);    
    // Limits of the integration
    etap = 9;    
    etam = FMAX(((*alfa) + sqrt(1 - (*tau) * (*tau)) * NIM8) / (*tau), -8);
    etamm = ((*alfa) + sqrt(1 - (*tau) * (*tau)) * (-NIM15)) / (*tau);
//    etam = FMAX(((*alfa) + sqrt(1. - (*tau) * (*tau)) * normsinv(1.e-7)) / (*tau), -7.);
//    etam = ((*alfa) + sqrt(1. - (*tau) * (*tau)) * normsinv(1.e-7)) / (*tau);

    if (etamm > etap)
    {
      status = gsl_integration_qng (&F, etam, etap, aerr * (etap - etam) / 2., 0., &f1, &error, &neval) ;
    }
    else
    {
      status = gsl_integration_qng (&F, etam, etamm, aerr * (etap - etam) / 2., 0., &f1, &error, &neval) ;
      f2 = SQRT2PI * (1 - erf(etamm)) / (*im); // Analythical
      f1 += f2;
    }
    
    return f1;
  }
}
