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

#ifndef __MYMATH_H
#define __MYMATH_H

#include <cmath>
#include <math.h>

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
        
#define  AA1  (-3.969683028665376e+01)
#define  AA2   2.209460984245205e+02
#define  AA3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  _D1   7.784695709041462e-03
#define  _D2   3.224671290700398e-01
#define  _D3   2.445134137142996e+00
#define  _D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575

namespace pmfpack
{
  typedef double (*simple_function)(double);  
  extern simple_function erfinv;
  extern simple_function erf;
  double normsinv(double);
  double Erf(double x);
  double Erfinv(double x);
  double erfinvC(double x);
  double eval_pol ( double a[], int *n, double *x );
  double stvaln ( double *p );
  double dinvnr ( double *p, double *q );
  void cdfnor ( int *which, double *p, double *q, double *x, double *mean,
  double *sd, int *status, double *bound );
  double erfC(double);
  double fifdint (double a);
  void cumnor (double *arg, double *result, double *ccum );  
}

#endif
