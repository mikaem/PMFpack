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

#include "PMFMath.h"

// The implemented error function and its inverse are defined using 
// the normal distribution and not the error function as defined in
// e.g. http://mathworld.wolfram.com/Erf.html. The correlation between 
// the regular error function and the implemented based on the normal 
// distribution is as folows
//
//  erf(z) = 0.5 * (1 + erf_(z/sqrt(2)))
//
//  erfinv(z) = sqrt(2) * erfinv_(2 * z - 1)
// 
//  where
//                              / x
//                     2       |       -t*t
//       erf_(x) = ----------- |      e       dt
//                 sqrt(pi)    |
//                             /0
//
//                              / x
//                     1       |       -t*t/2
//        erf(x) = ----------- |      e       dt
//                  sqrt(2 pi) |
//                             /-oo
// 
// and erfinv and erfinv_ are the inverses of erf and erf_ respectively

namespace pmfpack
{

  double gsl_erf(double z)
  {
    return 0.5 * (1 + gsl_sf_erf(z / M_SQRT2));
  }

  #ifdef HAS_BOOST
  // Boost use the erf_ definition and as such we need to wrap these up
  double boost_erf(double z)
  {
    return 0.5 * (1 + boost::math::erf(z / M_SQRT2));
  }

  boost::math::normal s;
  double boost_cdf(double z)
  {
    return boost::math::cdf(s, z);
  }

  double boost_erfinv(double z)
  {
    return M_SQRT2 * boost::math::erf_inv(2 * z - 1);
  }
  #endif

  // Choose provider of error functions. I usually find GSL to be fastest (M. M.)
    simple_function erfinv = &gsl_cdf_ugaussian_Pinv;
  //  simple_function erfinv = &boost_erfinv;
    simple_function erf = &gsl_cdf_ugaussian_P;
  //  simple_function erf = &boost_cdf;
  //  simple_function erf = &gsl_erf;
  //  simple_function erf = &boost_erf;
}