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
    double dh, tau_p, tau_m;

    // Find appropriate dh.    
    dh = gsl_min((*fmean), 1 - (*fmean));  //   0 < f < 1
    dh = gsl_min(dh, gsl_min((*sigma), (*fmean) * (1 - (*fmean)) - (*sigma)));            //   0 < sigma < f * (1 - f)
    dh = gsl_min(dh / 10., 1e-5);
    dh = pow(2, (int) log2(dh));
    
    // Compute central tau first using a robust bracketing algorithm 
    roots->froot->compute(0);
    
    tau_p = dfunction_df(+dh, roots);
    tau_m = dfunction_df(-dh, roots);
    (*dtaudf) = (tau_p - tau_m) / 2. / dh;
    (*d2taudfdf) = (tau_p - 2. * (*tau) + tau_m) / dh / dh;
    
    tau_p = dfunction_ds(+dh, roots);
    tau_m = dfunction_ds(-dh, roots);
    (*dtauds) = (tau_p - tau_m) / 2. / dh;
    (*d2taudsds) = (tau_p - 2. * (*tau) + tau_m) / dh / dh;
  }
}
