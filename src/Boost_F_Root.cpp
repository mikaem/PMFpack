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

#include "Boost_F_Root.h"

using namespace std;

namespace pmfpack{
    
  Boost_F_Root::Boost_F_Root(Integrator *_integrator)
  : BoostRoot(_integrator)
  {
    F = new BoostFunction(integrator, &func);
    tol = new tolerance(1.e-12);
    lower = 1.e-3;    
    upper = 1.-1.e-6;    
    fsolver = 0;
    maxiter = 100;
  }
    
  Boost_F_Root::~Boost_F_Root()
  {
    delete F;   
    delete tol;
  }

  void Boost_F_Root::realloc(const unsigned int f)
  {
    assert( 0 <= f <= 2);
    fsolver = f;
  }

  double Boost_F_Root::compute(int verbose)
  {
    static double x;
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
    if (fsolver == 0)
    {
      found = boost::math::tools::toms748_solve(*F, lower, upper, *tol, it);
      if (verbose) std::cout << " Boost toms748 " << std::endl;
      if (it == maxiter) std::cout << " Root not found Boost toms748 Fsolver" << std::endl;
    }
    else if (fsolver == 1)
    {
      found = boost::math::tools::bracket_and_solve_root(*F, sqrt(1 - 2 * (*tau)), 1.1, false, *tol, it);
      if (verbose) std::cout << " Boost bracket_and_solve_root " << std::endl;
      if (it == maxiter) std::cout << " Root not found Boost bracket_and_solve Fsolver" << std::endl;      
    }
    else if (fsolver == 2)
    {
      found = boost::math::tools::bisect(*F, lower, upper, *tol, it);
      if (verbose) std::cout << " Boost bisect " << std::endl;

    }
    x = (found.first + found.second) / 2.;
    (*tau) = (1 - x * x) / 2;
    if (verbose) std::cout << " Root = " << (*tau) << " iterations " << maxiter << std::endl;
    if (reversed)
    {
      (*fmean) = 1. - (*fmean);
      (*im) = (*fmean) * (*fmean) + (*sigma);
    }    
    return (*tau);
  }
}