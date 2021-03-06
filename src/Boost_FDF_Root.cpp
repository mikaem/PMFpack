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

#include "Boost_FDF_Root.h"

using namespace std;

namespace pmfpack{
    
  Boost_FDF_Root::Boost_FDF_Root(Integrator *_integrator)
  : BoostRoot(_integrator)
  {
    F2 = new BoostF2Function(integrator, &fdfunc);
    F3 = new BoostF3Function(integrator, &fd2func);
    lower = 1.e-3;    
    upper = 1.-1.e-6;    
    fdfsolver = 0;
    maxiter = 100;
//     digits = numeric_limits<double>::digits/2;
    digits = 48;
  }
    
  Boost_FDF_Root::~Boost_FDF_Root()
  {
    delete F2;
    delete F3;  
  }

  void Boost_FDF_Root::realloc(const unsigned int f)
  {
    assert( 0 <= f <= 2);
    fdfsolver = f;
  }
  
  void Boost_FDF_Root::set_digits(const uint d)
  {
    digits = d;
  }

  double Boost_FDF_Root::compute(int verbose)
  {
    double x;
    boost::uintmax_t it;
    bool reversed;    

    if ((*fmean) > 0.5)
    {
      reversed = true;
      (*fmean) = 1. - (*fmean);
      (*im) = (*fmean) * (*fmean) + (*sigma);
    }
    else{
      reversed = false;
    }

    (*alfa) = erfinv(1 - (*fmean));    
    it = maxiter;
    if (fdfsolver == 0)
    {
      x = boost::math::tools::newton_raphson_iterate(*F2, sqrt(1 - 2 * (*tau)), lower, upper, digits, it);
      if (verbose) cout << " Boost Newton Raphson " << endl;
    }
    else if (fdfsolver == 1)
    {
      x = boost::math::tools::halley_iterate(*F3, sqrt(1 - 2 * (*tau)), lower, upper, digits, it);
      if (verbose) cout << " Boost halley_iterate " << endl;
    }
    else if (fdfsolver == 2)
    {
      x = boost::math::tools::schroeder_iterate(*F3, sqrt(1 - 2 * (*tau)), lower, upper, digits, it);
      if (verbose) cout << " Boost schroeder_iterate " << endl;
    }
    (*tau) = (1 - x * x) / 2;
    if (verbose) cout << " Root = " << (*tau) << " iterations " << it << endl;
    if (reversed)
    {
      (*fmean) = 1. - (*fmean);
      (*im) = (*fmean) * (*fmean) + (*sigma);
      (*alfa) = erfinv(1 - (*fmean));
    }
    return (*tau);
  }
}