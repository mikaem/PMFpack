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

#include "Lookup.h"

using namespace std;

namespace pmfpack
{
  double ***f3tensor(int nrow,int ncol,int ndep)
  {
    double ***t;
    t = (double ***) malloc(nrow * sizeof(double**));
    t[0] = (double **) malloc((nrow * ncol) * sizeof(double*));
    t[0][0] = (double *) malloc((nrow * ncol * ndep) * sizeof(double));
    for(int j=1; j<ncol; j++) t[0][j] = t[0][j-1] + ndep;
    for(int i=1; i<nrow; i++) {
        t[i] = t[i-1] + ncol;
        t[i][0] = t[i-1][0] + ncol * ndep;
        for(int j=1; j<ncol; j++) t[i][j] = t[i][j-1] + ndep;
    }
    return t;
  }

  void free_f3tensor(double ***t)
  {
    free(t[0][0]);
    free(t[0]);
    free(t);
  }
  
  Lookup::Lookup(Derivator *_derivator)
  : derivator(_derivator)
  {
    fmean = derivator->fmean;
    sigma = derivator->sigma;
    alfa = derivator->alfa;
    tau  = derivator->tau;
    im  = derivator->im;
    dtau[0] = derivator->tau;
    dtau[1] = derivator->dtaudf;
    dtau[2] = derivator->dtauds;
    dtau[3] = derivator->d2taudfdf;
    dtau[4] = derivator->d2taudsds;
  }
    
  void Lookup::init(int N1, int N2)
  {
    Nf = N1;
    Ns = N2;
  }
  
  void Lookup::arange(double *y, double x0, double dx, int N)
  {
    for(int i = 0; i < N; i++)
        y[i] = x0 + dx * (i * 1.);
  }

  void Lookup::uniformgrid(double xmax, double ymax)
  {
    double dy, dx;
    
    dy = 2. * M_PI / (Ns + 1.);
    arange(fm, -M_PI + dy, dy, Ns);
    dx = 2. * M_PI / (Nf + 1.);
    arange(Is, -M_PI + dx, dx, Nf);
        
    for(int i = 0; i < Ns; i++)
      Is[i] = Is[i] * ymax;
    for(int i = 0; i < Nf; i++)
      fm[i] = fm[i] * xmax;
    
    grid_type = 0;
  }
    
  void Lookup::locate(uint *m, uint *n)
  {
    double ss = (*sigma) / (*fmean) / (1. - (*fmean));
    if (grid_type == 1)
    {
      *m = floor((Nf + 1) / 2. / M_PI * tan(atan2(M_PI, 1.) * (2. * (*fmean) - 1.)) + (Nf + 1) / 2.) -1;
      *n = floor((Ns + 1) / 2. / M_PI * tan(atan2(M_PI, 1.) * (2. * ss - 1.)) + (Ns + 1) / 2.) - 1;      
    }
    else
    {
      locate1(fm, Nf+1, (*fmean), m);
      locate1(Is, Ns+1, (*sigma) / (*fmean) / (1. - (*fmean)), n);
    }
    (*m)++; 
    (*n)++;
  }
  
  void Lookup::locate1(double xx[], uint n, double x, uint *j)
  {
    int ju,jm,jl,i;
    int ascnd;

    jl=0;
    ju=n;
    ascnd=(xx[n-2] >= xx[0]);
    while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
  //        jm = (ju+jl)/2;
        if ((x >= xx[jm-1]) == ascnd){
            jl=jm;}
        else{
            ju=jm;}   
    }
    if (x == xx[0]) *j=0;
    else if(x == xx[n-1]) *j=n-2;
    else *j=jl-1;
  }
  
  void Lookup::atangrid(double xmax, double ymax)
  {
    double dy, dx, at;
    
    at = atan2(M_PI, 1.);
    dy = 2. * M_PI / (Ns + 1.);
    arange(Is, -M_PI + dy, dy, Ns);
    for(int i = 0; i < Ns; i++)
      Is[i] = (atan2(Is[i], 1.) / at + 1.) / 2.;
    dx = 2. * M_PI / (Nf + 1.);
    arange(fm, -M_PI + dx, dx, Nf);
    for(int i = 0; i < Nf; i++)
      fm[i] = (atan2(fm[i], 1.) / at + 1.) / 2.;        
    for(int i = 0; i < Ns; i++)
      Is[i] = Is[i] * ymax;
    for(int i = 0; i < Nf; i++)
      fm[i] = fm[i] * xmax;
    
    grid_type = 1;
    
//     for(int i = 0; i < Nf; i++)
//         std::cout << fm[i] << " " << Is[i] << std::endl;
  }
  
}