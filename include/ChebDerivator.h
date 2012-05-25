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

#ifndef __CHEBDERIVATOR_H
#define __CHEBDERIVATOR_H

#include "Derivator.h"
#include <gsl/gsl_chebyshev.h>

namespace pmfpack
{

  class ChebDerivator : public Derivator 
  {
  public:
    
    ChebDerivator(Roots *, int);
            
    ~ChebDerivator();
      
    virtual double compute(int);
    
    int cheb_order;
    
    gsl_cheb_series *c, *c1, *c2;
    
    gsl_function F;

  };
}

#endif
