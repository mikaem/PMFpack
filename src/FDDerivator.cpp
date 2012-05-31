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

#include "FDDerivator.h"

namespace pmfpack
{
  FDDerivator::FDDerivator(Roots *_roots)
  : Derivator(_roots) {};
  
  double FDDerivator::compute(int verbose)
  {
    double fmean0, tau0, im0, dh, dh0;
    double tau_p, tau_m, s, s1,i1;
    int    dhi, count;

    fmean0 = (*fmean);
    im0 = (*im);

    // Set the s pointer to either sigma or im
    s = *roots->central ? (*sigma) : (*im); 
    // Get the largest possible dh. This is computed from 0<I<1, by taking all combinations of f+-dh,s+-dh.
    // This is valid also for d2tau/dfmean/dsigma
    if (!(*roots->central)){
      dh = sqrt(s) - fmean0;
      dh = gsl_min(dh, 0.5 * (fmean0 - s));
      dh = gsl_min(dh, s - fmean0 * fmean0);
      dh = gsl_min(dh, ((-(2. * fmean0 - 1.) + sqrt(4. * (s - fmean0) + 1.)) / 2.));
      dh = gsl_min(dh, ((2. * fmean0 + 1. + sqrt(4. * (s + fmean0) + 1.)) / 2.));
      dh = gsl_min(dh, ((-(2. * fmean0 + 1.) + sqrt(4. * (s + fmean0) + 1.)) / 2.));   
      dh = gsl_min(dh, ((2. * fmean0 - 1. + sqrt(4. * (s - fmean0) + 1.)) / 2.));
    }
    else{
      dh = (-(2. * fmean0 - 1.) + sqrt(1. - 4. * s)) / 2.;
      dh = gsl_min(dh, (2. * fmean0 - 1. + sqrt(1. - 4. * s)) / 2);
      dh = gsl_min(dh, fmean0 - (fmean0 * fmean0 + s));
      dh = gsl_min(dh, -fmean0 + sqrt(fmean0 - s));
      dh = gsl_min(dh, fmean0 - 1. + sqrt(1. - fmean0 - s));
      dh = gsl_min(dh, -(fmean0 - 1.) + sqrt(1. - fmean0 - s));
      dh = gsl_min(dh, s);
    }
    dh = gsl_min(5.e-4, dh / 250.);
    dhi = (int) log10(dh);
    dh = pow(10, dhi);
    dh0 = dh;
    
    // Compute central tau first using a robust bracketing algorithm 
    roots->froot->compute(0);
    tau0 = (*tau);
    tau_m = 1.;
    count = 0;
    
    while((tau_m > 0.005 || tau_m < 0.0005) && count < 5){
      // Find an appropriate dh
      if (tau_m > 0.005) dh = dh / 10;
      else dh = dh * 10;
      tau_p = dfunction_df(+dh, roots);
      tau_m = (fabs(tau_p - tau0) / tau0);
      count++;
    }
    tau_m = dfunction_df(-dh, roots);
    (*dtaudf) = (tau_p - tau_m) / 2. / dh;
    (*d2taudfdf) = (tau_p - 2. * tau0 + tau_m) / dh / dh;
    
    count = 0;
    while((tau_m > 0.005 || tau_m < 0.0005) && count < 5){
      // Find an appropriate dh
      if (tau_m > 0.005) dh = dh / 10;
      else dh = dh * 10;
      tau_p = dfunction_ds(+dh, roots);
      tau_m = (fabs(tau_p - tau0) / tau0);
      count++;
    }
    tau_m = dfunction_ds(-dh, roots);
    (*dtauds) = (tau_p - tau_m) / 2. / dh;
    (*d2taudsds) = (tau_p - 2. * tau0 + tau_m) / dh / dh;
    
  }
}
