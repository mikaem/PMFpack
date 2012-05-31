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
// First added:  2012-05-23
// Last changed: 2012-05-23

#ifndef __GSLLOOKUP_H
#define __GSLLOOKUP_H

#include "Lookup.h"
#include "gsl/gsl_interp.h"
#include <stdio.h>
#include <stdlib.h>

namespace pmfpack
{
  
  class GSLLookup : public Lookup
  {
  public:
                      
    GSLLookup(Derivator *);
    
    ~GSLLookup();    
                
    void compute(int, bool derivatives=false, bool polish=false);
    
    void operator()(int tau, bool derivatives=false, bool polish=false) 
    {
      compute(tau, derivatives, polish);
    }
    
    void init(int, int);
    
    int read_table(char*);
    
    void generate_table(int, int, char*);
    
    void set_order(int);
    
    void lookup_clean();
    
    void gslpolin2(double [], double [], double **, double *);
    
  private:
    // Pointers to the data in lookuptables. 
    double ***tau_table;
    
    double *ys[16];  // Submatrix used by interpolation routine
    
    const gsl_interp_type *t;
    
    gsl_interp *s;
    
    int order;  // order of interpolation

  };

}

#endif
