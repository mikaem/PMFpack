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

namespace pmfpack
{
  
  double gslfunction(double x, void *params)
  {
    double **data = (double **) params;
    return function(x, data);    
  }

  GSLIntegrator::GSLIntegrator(double **_parameters)
  : Integrator(_parameters)
  {
    F.function = &gslfunction;    
    F.params = parameters;
    w = gsl_integration_workspace_alloc (1000);      
    table = gsl_integration_glfixed_table_alloc(96);   
  }
  
  GSLIntegrator::~GSLIntegrator()
  {
    gsl_integration_workspace_free(w);
    gsl_integration_glfixed_table_free(table);
  }
  
  double GSLIntegrator::compute()
  {
    int status;
    double phi_min, phi_X_unity, phi_max, integral, aerr, error;    
    size_t neval;
    
    // Note. We are solving for x = sqrt(1 - 2 * tau). tau below is thus really x and not tau as given in the papers
    
    aerr = 5.e-8 / (SQRT2PI);    
    // Limits of the integration
    phi_max = 9;    
    phi_min = FMAX(((*alfa) + sqrt(1 - (*tau) * (*tau)) * NIM8) / (*tau), -9);
    phi_X_unity = ((*alfa) + sqrt(1 - (*tau) * (*tau)) * (-NIM15)) / (*tau); // above this X in integrand is > 1-1e-15 and can be set to unity, which makes the integral analytical

    if (phi_X_unity > phi_max)
    {
      //status = gsl_integration_qng (&F, phi_min, phi_max, aerr * (phi_max - phi_min) / 2., 0., &integral, &error, &neval) ;
      //status = gsl_integration_qag (&F, phi_min, phi_max, aerr * (phi_max - phi_min) / 2., 0., w->limit, 2, w, &integral, &error) ;
      integral = gsl_integration_glfixed (&F, phi_min, phi_max, table);
    }
    else
    {
      //status = gsl_integration_qng (&F, phi_min, phi_X_unity, aerr * (phi_max - phi_min) / 2., 0., &integral, &error, &neval) ;
      integral = gsl_integration_glfixed (&F, phi_min, phi_max, table);
      integral += SQRT2PI * (1 - erf(phi_X_unity)); // Analytical
    }
    
    return integral / (*im);
  }
}
