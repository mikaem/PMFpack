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

#include "Integrator.h"
#include "GSLIntegrator.h"
#include <iostream>
#include <iomanip>

namespace pmfpack
{
  double func(double, Integrator*);
  double dfunc(double, Integrator*);
  void fdfunc(double, Integrator*, double *, double *);
  void fd2func(double, Integrator*, double *, double *, double *);
  
  class Root
  {
  public:

    Root(Integrator *);
    
    virtual ~Root() {};    
    
    virtual double compute(int) = 0;
    
    virtual void realloc(const unsigned int) = 0;
    
    virtual void error_message_off() = 0;
            
    int max_iteration;
    
    double **parameters;
    
    double *fmean, *sigma, *alfa, *tau, *im, *f, *df;
    
    bool central;
    
    Integrator *integrator;
      
  };
  
  // Create a struct that holds both an froot and an fdfroot class
  typedef struct _Roots
  {
    Root *froot, *fdfroot;
    bool *central;
  } Roots;

}

#endif
