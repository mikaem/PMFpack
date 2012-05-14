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

#ifndef __ROOT_H
#define __ROOT_H

#include <iostream>
#include <iomanip>
#include "Integrator.h"
#include "GSLIntegrator.h"
#include "MyMath.h"

namespace pmfpack
{
  double func(double, void *);
  
  double dfunc(double, void *);
  
  void fdfunc(double, void *, double *, double *);
  
  class Root
  {
  public:
                      
    Root();
    
    Root(Integrator *);
    
    ~Root() {};    
    
    virtual double compute(int) = 0;
    
    virtual void realloc(const unsigned int) = 0;
    
    virtual void error_message_off() = 0;
            
    int max_iteration;
    
    double **parameters;
    
    double *fmean, *sigma, *alfa, *tau, *im, *f, *df;
    
    bool central;
    
    Integrator *integrator;
      
  };

}

#endif
