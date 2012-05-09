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

#ifndef __GSLINTEGRATOR_H
#define __GSLINTEGRATOR_H

#include "Integrator.h"
#include <gsl/gsl_integration.h>
#include "gsl/gsl_cdf.h"

namespace pmfpack
{

  class GSLIntegrator : public Integrator 
  {
  public:
      
    GSLIntegrator(double **);
    
    ~GSLIntegrator();
    
    gsl_function F;
        
    gsl_integration_workspace *w;
    
    virtual double compute(const double*, const double*, const double*);
                      
  };
}

#endif
