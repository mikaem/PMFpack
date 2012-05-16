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
    F = new BoostF2Function(integrator, &fdfunc);
    lower = 1.e-3;    
    upper = 1.-1.e-6;    
    fdfsolver = 0;
    maxiter = 100;
    digits = std::numeric_limits<double>::digits;
  }
    
  Boost_FDF_Root::~Boost_FDF_Root()
  {
    delete F;  
  }

  void Boost_FDF_Root::realloc(const unsigned int f)
  {
    assert( 0 <= f <= 2);
    fdfsolver = f;
  }

  double Boost_FDF_Root::compute(int verbose)
  {
    double x;

    (*alfa) = erfinv(1 - (*fmean));    
    if (fdfsolver == 0)
    {
      x = boost::math::tools::newton_raphson_iterate(*F, sqrt(1 - 2 * (*tau)), lower, upper, digits);
      if (verbose) std::cout << " Boost Newton Raphson " << std::endl;
    }
    else if (fdfsolver == 1)
    {
//       found = boost::math::tools::bracket_and_solve_root(*F, sqrt(1 - 2 * (*tau)), 1.1, false, *tol, maxiter);
//       if (verbose) std::cout << " Boost bracket_and_solve_root " << std::endl;
    }
    else if (fdfsolver == 2)
    {
//       found = boost::math::tools::bisect(*F, lower, upper, *tol, maxiter);
//       if (verbose) std::cout << " Boost bisect " << std::endl;

    }
    (*tau) = (1 - x * x) / 2;
    if (verbose) std::cout << " Root = " << (*tau) << " iterations " << maxiter << std::endl;
    
    return (*tau);
  }
}