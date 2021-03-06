// Copyright (C) 2012 Mikael Mortensen
//
// This file is part of PMFpack.
//
// PMFpack is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNs later version.
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

#include "GSLLookup.h"

namespace pmfpack
{
  
  GSLLookup::GSLLookup(Derivator *_derivator)
  : Lookup(_derivator) 
  {
    s = NULL;    
    fm = NULL;
    t = gsl_interp_polynomial;
    set_order(4);
  }  
  
  void GSLLookup::init(int N1, int N2)
  {
    Lookup::init(N1, N2);
    fm = (double *) malloc(Nf * sizeof(double));
    Is = (double *) malloc(Ns * sizeof(double));
    atangrid(1.0, 1.0);
    tau_table = f3tensor(5, Nf, Ns);
  } 
  
  GSLLookup::~GSLLookup()
  {
    lookup_clean();
  }
  
  void GSLLookup::lookup_clean()
  {
    free(fm);
    free(Is);
    free_f3tensor(tau_table);
    gsl_interp_free(s);
  }
  
  int GSLLookup::read_table(char* filename)
  {
    FILE *stream;
    int status;
        
    stream = fopen(filename, "r");
    if ((stream = fopen(filename, "r")) == NULL)
      return -1;
    status = fscanf(stream, "%d", &Nf);
    status = fscanf(stream, "%d", &Ns);

    if (!fm) init(Nf, Ns);
    for (int i = 0; i < Nf; i++)
      status = fscanf(stream, "%lf", &fm[i]);
    for (int i = 0; i < Ns; i++)
      status = fscanf(stream, "%lf", &Is[i]);
    status = fscanf(stream, "%d", &grid_type);
        
    for (int k = 0; k < 5; k++){   
      for (int i = 0; i < Nf; i++){
        for (int j = 0; j < Ns; j++){  
          status = fscanf(stream, "%lf", &tau_table[k][i][j]);
        }
      }
    }    
    fclose(stream);
    
//     std::cout << " Read table" << std::endl;    
//     for (uint i=0; i<5; i++)
//     {
//       for (uint j=0; j<Nf; j++)
//       {
//         for (uint k=0; k<Ns; k++)
//           std::cout << i << " " << j << " " << k << " " << tau_table[i][j][k] << std::endl;
//       }
//     }
// 
    
    return 1;
  }
  
  void GSLLookup::set_order(int _order)
  {
    if (s != NULL)
      gsl_interp_free(s);
    order = _order;
    s = gsl_interp_alloc(t, order);
  }
  
  void GSLLookup::generate_table(int N1, int N2, char *filename)
  {
    double fmean0, sigma0, alfa0, tau0, im0;
    double err;
    FILE *f1;
    
    if (fm != NULL && (N1 != Nf || N2 != Ns))
    {
        lookup_clean();
        init(N1, N2);
    }        
    else if (fm == NULL) 
    {
        init(N1, N2);
    }
    // Remember original data
    fmean0 = (*fmean);
    sigma0 = (*sigma);
    alfa0  = (*alfa);
    tau0   = (*tau);
    im0    = (*im);
    
    // Open lookuptables for writing
    if ((f1 = fopen(filename, "w")) == NULL)
      fprintf(stderr, "Data file %s could not be opened\n", filename);
    
    fprintf(f1, "%d %d\n", Nf, Ns);
    for(int i = 0; i < Nf; i++)
      fprintf(f1, "%4.16e\n", fm[i]);
    for(int i = 0; i < Ns; i++)
      fprintf(f1, "%4.16e\n", Is[i]);
    fprintf(f1, "%d\n", grid_type);
    for (int k = 0; k < 5; k++){
      for(int i = 0; i < Nf; i++){
        for(int j = 0; j < Ns; j++){
          (*fmean) = fm[i];
          (*sigma) = Is[j] * fm[i] * (1. - fm[i]);
          (*im) = (*sigma) + fm[i] * fm[i];
          derivator->compute(0);        
          tau_table[k][i][j] = (*dtau[k]);
          fprintf(f1,"%4.16e\n", (*dtau[k]));
        }
      }
    }
    fclose(f1);
    
    // Reset data
    (*fmean) = fmean0;
    (*sigma) = sigma0;
    (*alfa)  = alfa0;
    (*tau)   = tau0;
    (*im)    = im0;
  }   

  void GSLLookup::gslpolin2(double x1a[], double x2a[], double **ya, double *y)
  {
    double *ymtmp = (double*) malloc(order * sizeof(double));
    for (int j = 0; j < order; j++) { 
      gsl_interp_init(s, x2a, ya[j], order);
      ymtmp[j] = gsl_interp_eval(s, x2a, ya[j], *sigma / *fmean / (1. - (*fmean)), NULL);
    }
    gsl_interp_init(s, x1a, ymtmp, order);
    *y = gsl_interp_eval(s, x1a, ymtmp, *fmean, NULL);
    free(ymtmp);
  }
  
  void GSLLookup::compute(int verbose, bool derivatives, bool polish)
  {
    uint m, n, ms1, ms2, ns1, ns2, ms, ns, check;
    double *x1s, *x2s, dy;

    // Locate point in lookuptable
    locate(&m, &n);
    
    // Compute if point is outside of table borders
    if (m == 0 || n == 0 || m == Nf || n == Ns)
    {
      if (derivatives)
        derivator->compute(verbose);
      else
        derivator->roots->froot->compute(verbose);
      return;
    }
    // Locate submatrix
    // Linear interpolation requires at least two points
    ms1 = IMIN(m+2, order), ms2 = IMAX(2, Nf-m+2), (ms1 < ms2)  ?  (ms=ms1) : (ms=ms2);
    ns1 = IMIN(n+2, order), ns2 = IMAX(2, Ns-n+2), (ns1 < ns2)  ?  (ns=ns1) : (ns=ns2);
    
    if (ms1 < order){
      ms2 = 2 * order - ms1;
    }
    else if(ms2 < 4){
      ms1 = 2 * order - ms2;
    }
    else{
      ms1 = order;
      ms2 = order;
    }
    if(ns1 < order){
      ns2 = 2 * order - ns1;
    }
    else if(ns2 < 4){
      ns1 = 2 * order - ns2;
    }
    else{
      ns1 = order;
      ns2 = order;
    }
    ms1 = (int)(floor(ms1 / 2.));
    ms2 = (int)(ceil(ms2 / 2.));
    ns1 = (int)(floor(ns1 / 2.));
    ns2 = (int)(ceil(ns2 / 2.));
    x1s = &fm[m-ms1];
    x2s = &Is[n-ns1];
    
    // Copy submatrix to ys and interpolate
    for(int i = m-ms1; i < m+ms2; i++){
      ys[i-m+ms1] = tau_table[0][i]+n-ns1;
      std::cout << *ys[i-m+ms1] << " " << i-m+ms1 << " " << i << " " << n-ns1 << std::endl;      
    }
    gslpolin2(x1s, x2s, ys, tau);
    // Polish the root with an fdf solver for best possible accuracy. The lookup is then initial guess.
    if (polish){
      if (derivatives)
        derivator->compute(verbose);
      else
        derivator->roots->fdfroot->compute(verbose);
      return; // Return because you don't need to look up derivatives you have already computed
    }
    if (derivatives){
      for (int k = 1; k < 5; k++){
        for(int i = m-ms1; i < m+ms2; i++){
          ys [i-m+ms1] = tau_table[k][i]+n-ns1;
        }
        gslpolin2(x1s, x2s, ys, dtau[k]);
      }    
    }
  }

  double ** GSLLookup::get_table(int i)
  {
    return tau_table[i];
  }
}
