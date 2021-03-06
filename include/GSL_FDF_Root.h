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

#ifndef __GSL_FDF_ROOT_H
#define __GSL_FDF_ROOT_H

#include "GSLRoot.h"

namespace pmfpack
{
  gsl_root_fdfsolver * gsl_root_fdfsolver_realloc (gsl_root_fdfsolver * s, const gsl_root_fdfsolver_type * T);

  class GSL_FDF_Root : public GSLRoot
  {
  public:
    GSL_FDF_Root(Integrator *);
    
    ~GSL_FDF_Root();
            
    gsl_function_fdf F;
    
    double compute(int);
    
    void realloc(const unsigned int);
    
    const gsl_root_fdfsolver_type *T;
    
    gsl_root_fdfsolver *s;
    
  };  

}

#endif