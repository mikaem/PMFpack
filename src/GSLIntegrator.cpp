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
#include "MyMath.h"
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
  
  double GSLIntegrator::compute(const double* etam,const double* etamm,const double* etap)
  {
    int status;
    size_t neval;
    double f1,f2,sr,aerr,error;    
    aerr = 5.e-8 / (SQRT2PI);

    if((*etamm) > (*etap))
    {
      status = gsl_integration_qng (&F,*etam,*etap, aerr*((*etap)-(*etam))/2.,0.,&f1, &error, &neval) ;
    }
    else
    {
      status = gsl_integration_qng (&F, *etam, *etamm, aerr*((*etap) - (*etam)) / 2., 0., &f1, &error, &neval) ;
      f2 = SQRT2PI * (1 - gsl_cdf_ugaussian_P(*etamm)) / (*parameters[4]); // Analythical
      f1 = f1 + f2;
    }
    
    return f1;
  }
}
