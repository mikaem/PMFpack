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

#ifndef __LOOKUP_H
#define __LOOKUP_H

#include "Derivator.h"

namespace pmfpack
{
  double ***f3tensor(int, int, int);  
  void free_f3tensor(double ***);
  
  class Lookup
  {
  public:
                      
    Lookup(Derivator *);
    
    virtual ~Lookup() {};    
        
    Derivator *derivator;
    
    double *fmean, *sigma, *alfa, *tau, *im, *f, *df;
    
    double* dtau[5]; // tau, dtaudf, dtauds, d2taudfdf, d2taudsds
    
    int Nf, Ns; // Number of points in fmean and sigma space
    
    double *fm, *Is;
    
    uint grid_type;
    
    void uniformgrid(double, double);
    
    void atangrid(double, double);
    
    void arange(double *, double, double, int);
    
    void locate1(double [], uint, double, uint *);
    
    void locate(uint *, uint *);
        
    virtual void init(int, int);
    
    virtual void operator()(int, bool derivatives=false, bool polish=false) = 0;
    
    virtual void compute(int, bool derivatives=false, bool polish=false) = 0;
    
    virtual int read_table(char*) = 0;
    
    virtual void generate_table(int, int, char*) = 0;
    
    virtual double ** get_table(int) = 0;
    
  };

}

#endif
