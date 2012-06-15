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

#ifndef __PMFMATH_H
#define __PMFMATH_H

#ifdef HAS_BOOST
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#endif
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
        
namespace pmfpack
{
  static double sqrarg;
  #define PI 3.1415926535897931018
  #define SQRT2 1.4142135623730952219
  #define SQRTPI 1.77245385090551590643
  #define SQRT2PI 2.50662827463100020000
  #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
  #define SQRT(a) ((sqrarg=(a)) == 0.0 ? 0.0 : pow(sqrarg,0.5))
  // Inverse error functions NIM7=Erfinv(1e-7)
  #define NIM7 -5.1993375821928174
  #define NIM8 -5.6120012441747891
  #define NIM9 -5.9978070150076865
  #define NIM14 -7.6506280929352686
  #define NIM15 -7.941345326170997
  #define NIM16 -8.2220822161304365
  static double maxarg1,maxarg2;
  #define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
          (maxarg1) : (maxarg2))
  #define FMIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) < (maxarg2) ?\
          (maxarg1) : (maxarg2))
  #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

  static unsigned int maxa1,maxa2;
  #define IMAX(a,b) (maxa1=(a),maxa2=(b),(maxa1) > (maxa2) ?\
          (maxa1) : (maxa2))
  #define IMIN(a,b) (maxa1=(a),maxa2=(b),(maxa1) < (maxa2) ?\
          (maxa1) : (maxa2))

  double gsl_erf(double);
  double boost_erf(double);
  double boost_cdf(double);
  double boost_erfinv(double);
  typedef double (*simple_function)(double);  
  extern simple_function erfinv;
  extern simple_function erf;
}

#endif
