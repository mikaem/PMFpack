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

// There is a slight difference between using central or integer moments,
// because only the integer is involved in the computation (see function 
// in Integrator.cpp).
// If we are using central moments in the CFD calculation we change fmean 
// in the finite difference approximations in this routine and keep the 
// central moment constant. Hence the integer that is used in computations 
// must be modified.
// Using integer moments in the CFD there's no need to do anything special
// with the finite difference calculations, since the integer moments 
// should be kept constant when altering fmean.

#include "Derivator.h"

namespace pmfpack
{
  Derivator::Derivator(bool _central, Root *_froot, Root *_fdfroot)
  {
    params = new Params();
    params->central = _central;
    params->froot = _froot;
    params->fdfroot = _fdfroot;
    fmean = _froot->parameters[0];
    sigma = _froot->parameters[1];
    alfa  = _froot->parameters[2];
    tau   = _froot->parameters[3];
    im    = _froot->parameters[4];
    dtaudf = _froot->parameters[7];
    dtauds = _froot->parameters[8];
    d2taudfdf = _froot->parameters[9];
    d2taudsds = _froot->parameters[10];
  }
  
  double dfunction_df(const double df, void *pars)
  {
    /*
     * Add a tiny df to fmean, modify im, compute tau and then reset fmean + im.
     */
    double fmean0, im0, alfa0, tau0, newtau;
    Params *params = (Params *) pars;
    double **p = params->fdfroot->parameters;
    
    // Remember old values
    fmean0 = (*p[0]);
    alfa0 = (*p[2]);
    tau0 = (*p[3]);
    im0 = (*p[4]);
        
    // Modify parameters
    (*p[0]) = (*p[0]) + df;
    (*p[4]) = params->central ? (*p[0]) * (*p[0]) + (*p[1]) : (*p[4]);
    
    // Compute new tau
    params->fdfroot->compute(0);
    
    // Reset to old values
    (*p[0]) = fmean0;
    (*p[2]) = alfa0;
    (*p[4]) = im0;
    newtau = (*p[3]);
    (*p[3]) = tau0;

    return newtau;
  }
  
  double dfunction_ds(const double ds, void *pars)
  {
    /*
     * Add a tiny ds to variance, modify im, compute tau and then reset im.
     */
    double im0, tau0, newtau;
    Params *params = (Params *) pars;
    double **p = params->fdfroot->parameters;
    
    // Remember old values
    tau0 = (*p[3]);
    im0 = (*p[4]);
        
    // Modify parameters (only im used in fdfroot->compute)
    (*p[4]) = (*p[4]) + ds;

    // Compute new tau
    params->fdfroot->compute(0);
    
    // Reset to old values
    (*p[4]) = im0;
    newtau = (*p[3]);
    (*p[3]) = tau0;

    return newtau;
  }
}
