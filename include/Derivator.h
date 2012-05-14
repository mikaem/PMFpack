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
// First added:  2012-05-11
// Last changed: 2012-05-11

#ifndef __DERIVATOR_H
#define __DERIVATOR_H

#include "Root.h"

namespace pmfpack
{
  struct Params
  {
    Root *froot, *fdfroot;
    bool central;
  };
  
  double dfunction_df(const double, void *);
  double dfunction_ds(const double, void *);

  class Derivator
  {
  public:
    Derivator(bool, Root *, Root *);
            
    ~Derivator() {};        
    
    virtual double compute(int) = 0;
    
    double *fmean, *sigma, *alfa, *tau, *im;
    double *dtaudf, *dtauds, *d2taudfdf, *d2taudsds;
    
    Params *params;
        
  };
}

#endif
