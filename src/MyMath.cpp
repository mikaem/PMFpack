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

#include "MyMath.h"
#include <iostream>

namespace pmfpack
{
      
  double normsinv(double p)
  {
    double x;
    double q,r,e,u;
    
//         std::cout << "  " << p << " " << P_LOW << std::endl; 
    if ((0 < p )  && (p < P_LOW))
    {
      q = sqrt(-2*log(p));
      x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((_D1*q+_D2)*q+_D3)*q+_D4)*q+1);
    }
    else
    {
      if ((P_LOW <= p) && (p <= P_HIGH))
      {
        q = p - 0.5;
        r = q*q;
        x = (((((AA1*r+AA2)*r+AA3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
      }
      else{
        if ((P_HIGH < p) && (p < 1))
        {
          q = sqrt(-2*log(1-p));
          x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((_D1*q+_D2)*q+_D3)*q+_D4)*q+1);
        }
      }
    }

    // Unblock for machine precision

    //if(( 0 < p)&&(p < 1)){
    //   e = Erf(x) - p;
    //   u = e * SQRT2PI * exp(x*x/2);
    //   x = x - u/(1 + x*u/2);
    //}

    //

    return x;
  }
  
  void cumnor ( double *arg, double *result, double *ccum )

  /****************************************************************************
  // 
  //  Purpose:
  // 
  //    CUMNOR computes the cumulative normal distribution.
  //
  //  Discussion:
  //
  //    This function evaluates the normal distribution function:
  //
  //                              / x
  //                     1       |       -t*t/2
  //          P(x) = ----------- |      e       dt
  //                 sqrt(2 pi)  |
  //                             /-oo
  //
  //    This transportable program uses rational functions that
  //    theoretically approximate the normal distribution function to
  //    at least 18 significant decimal digits.  The accuracy achieved
  //    depends on the arithmetic system, the compiler, the intrinsic
  //    functions, and proper selection of the machine-dependent
  //    constants.
  //
  //  Author: 
  //
  //    W J Cody
  //    Mathematics and Computer Science Division
  //    Argonne National Laboratory
  //    Argonne, IL 60439
  //
  //  Reference:
  //
  //    W J Cody,
  //    Rational Chebyshev approximations for the error function,
  //    Mathematics of Computation, 
  //    1969, pages 631-637.
  //
  //    W J Cody, 
  //    Algorithm 715: 
  //    SPECFUN - A Portable FORTRAN Package of Special Function Routines 
  //      and Test Drivers,
  //    ACM Transactions on Mathematical Software,
  //    Volume 19, 1993, pages 22-32.
  //
  //  Parameters:
  //
  //    Input, double *ARG, the upper limit of integration.
  //
  //    Output, double *CUM, *CCUM, the Normal density CDF and
  //    complementary CDF.
  //
  //  Local Parameters:
  //
  //    Local, double EPS, the argument below which anorm(x) 
  //    may be represented by 0.5D+00 and above which  x*x  will not underflow.
  //    A conservative value is the largest machine number X
  //    such that   1.0D+00 + X = 1.0D+00   to machine precision.
  /*/
  {
  static double a[5] = {
      2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
      1.8154981253343561249e04,6.5682337918207449113e-2
  };
  static double b[4] = {
      4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
      4.5507789335026729956e04
  };
  static double c[9] = {
      3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
      5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
      1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
  };
  static double d[8] = {
      2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
      6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
      3.8912003286093271411e04,1.9685429676859990727e04
  };
  static double half = 0.5e0;
  static double p[6] = {
      2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
      1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
  };
  static double one = 1.0e0;
  static double q[5] = {
      1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
      3.78239633202758244e-3,7.29751555083966205e-5
  };
  static double sixten = 1.60e0;
  static double sqrpi = 3.9894228040143267794e-1;
  static double thrsh = 0.66291e0;
  static double root32 = 5.656854248e0;
  static double zero = 0.0e0;
  static int K1 = 1;
  static int K2 = 2;
  static int i;
  double del,eps,temp,x,xden,xnum,y,xsq,min;
  //
  //  Machine dependent constants
  //
  //    eps = dpmpar(&K1)*0.5e0;
  //    min = dpmpar(&K2);
      eps = 1.e-16;
      min = 0.5e-16;
      x = *arg;
      y = fabs(x);

      if(y <= thrsh) {
  //
  //  Evaluate  anorm  for  |X| <= 0.66291
  //
          xsq = zero;
          if(y > eps) xsq = x*x;
          xnum = a[4]*xsq;
          xden = xsq;
          for ( i = 0; i < 3; i++ )
          {
              xnum = (xnum+a[i])*xsq;
              xden = (xden+b[i])*xsq;
          }
          *result = x*(xnum+a[3])/(xden+b[3]);
          temp = *result;
          *result = half+temp;
          *ccum = half-temp;
      }
  //
  //  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
  //
      else if(y <= root32) {
          xnum = c[8]*y;
          xden = y;
          for ( i = 0; i < 7; i++ )
          {
              xnum = (xnum+c[i])*y;
              xden = (xden+d[i])*y;
          }
          *result = (xnum+c[7])/(xden+d[7]);
          xsq = fifdint(y*sixten)/sixten;
          del = (y-xsq)*(y+xsq);
          *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
          *ccum = one-*result;
          if(x > zero) {
              temp = *result;
              *result = *ccum;
              *ccum = temp;
          }
      }
  //
  //  Evaluate  anorm  for |X| > sqrt(32)
  //
      else  {
          *result = zero;
          xsq = one/(x*x);
          xnum = p[5]*xsq;
          xden = xsq;
          for ( i = 0; i < 4; i++ )
          {
              xnum = (xnum+p[i])*xsq;
              xden = (xden+q[i])*xsq;
          }
          *result = xsq*(xnum+p[4])/(xden+q[4]);
          *result = (sqrpi-*result)/y;
          xsq = fifdint(x*sixten)/sixten;
          del = (x-xsq)*(x+xsq);
          *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
          *ccum = one-*result;
          if(x > zero) {
              temp = *result;
              *result = *ccum;
              *ccum = temp;
          }
      }
      if(*result < min) *result = 0.0e0;
  //
  //  Fix up for negative argument, erf, etc.
  //
      if(*ccum < min) *ccum = 0.0e0;
      
  }
  //****************************************************************************
  double fifdint ( double a )

  //****************************************************************************
  //
  //  Purpose:
  // 
  //    FIFDINT truncates a double number to an integer.
  //
  //  Parameters:
  //
  // a     -     number to be truncated 
  {
  return (double) ((int) a);
  }

  double erfC(double x)
  {
    double p,q,zz;
    zz=SQRT2*x;
    cumnor(&zz,&p,&q);
    return 2.*p-1.;
  }

  double Erf(double x)
  {
    double p,q;
    int isn;
    isn=isinf(x);
    if(isn==1)return 1.;
    else if(isn==-1)return 0.;
    else{
    cumnor(&x,&p,&q);
    return p;}
  }

  double Erfinv(double x)
  {
    double pn,qn,xn,bound,sd,mean;
    int which,status;
    pn=x;
    qn=1.-pn;
    xn=-1;
    sd=1.;
    mean=0.;
    which=2;
    cdfnor ( &which, &pn, &qn, &xn, &mean, &sd, &status, &bound );
    return xn;
  }

  double erfinvC(double x)
  {
    double pn,qn,xn,bound,sd,mean;
    int which,status;
    pn=0.5*(1.+x);
    qn=1.-pn;
    xn=-1;
    sd=1.;
    mean=0.;
    which=2;
    cdfnor ( &which, &pn, &qn, &xn, &mean, &sd, &status, &bound );
    return xn/SQRT2;
  }

  void cdfnor ( int *which, double *p, double *q, double *x, double *mean,
  double *sd, int *status, double *bound )

  //****************************************************************************
  //
  //  Purpose:
  // 
  //    CDFNOR evaluates the CDF of the Normal distribution.
  //
  //  Discussion:
  //
  //    A slightly modified version of ANORM from SPECFUN
  //    is used to calculate the cumulative standard normal distribution.
  //
  //    The rational functions from pages 90-95 of Kennedy and Gentle
  //    are used as starting values to Newton's Iterations which 
  //    compute the inverse standard normal.  Therefore no searches are
  //    necessary for any parameter.
  //
  //    For X < -15, the asymptotic expansion for the normal is used  as
  //    the starting value in finding the inverse standard normal.
  //
  //    The normal density is proportional to
  //    exp( - 0.5D+00 * (( X - MEAN)/SD)**2)
  //
  //  Reference:
  //
  //    Abramowitz and Stegun,  
  //    Handbook of Mathematical Functions 
  //    1966, Formula 26.2.12.
  //
  //    W J Cody,
  //    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
  //      Special Function Routines and Test Drivers,
  //    ACM Transactions on Mathematical Software,
  //    Volume 19, pages 22-32, 1993.
  //
  //    Kennedy and Gentle,
  //    Statistical Computing,
  //    Marcel Dekker, NY, 1980,
  //    QAA276.4  K46
  //
  //  Parameters:
  //
  //    Input, int *WHICH, indicates which argument is to be calculated
  //    from the others.
  //    1: Calculate P and Q from X, MEAN and SD;
  //    2: Calculate X from P, Q, MEAN and SD;
  //    3: Calculate MEAN from P, Q, X and SD;
  //    4: Calculate SD from P, Q, X and MEAN.
  //
  //    Input/output, double *P, the integral from -infinity to X 
  //    of the Normal density.  If this is an input or output value, it will
  //    lie in the range [0,1].
  //
  //    Input/output, double *Q, equal to 1-P.  If Q is an input
  //    value, it should lie in the range [0,1].  If Q is an output value,
  //    it will lie in the range [0,1].
  //
  //    Input/output, double *X, the upper limit of integration of 
  //    the Normal density.
  //
  //    Input/output, double *MEAN, the mean of the Normal density.
  //
  //    Input/output, double *SD, the standard deviation of the 
  //    Normal density.  If this is an input value, it should lie in the
  //    range (0,+infinity).
  //
  //    Output, int *STATUS, the status of the calculation.
  //    0, if calculation completed correctly;
  //    -I, if input parameter number I is out of range;
  //    1, if answer appears to be lower than lowest search bound;
  //    2, if answer appears to be higher than greatest search bound;
  //    3, if P + Q /= 1.
  //
  //    Output, double *BOUND, is only defined if STATUS is nonzero.
  //    If STATUS is negative, then this is the value exceeded by parameter I.
  //    if STATUS is 1 or 2, this is the search bound that was exceeded.
  //
  {
  static int K1 = 1;
  static double z,pq;


  //
  //     Check arguments
  //
      *status = 0;
      if(!(*which < 1 || *which > 4)) goto S30;
      if(!(*which < 1)) goto S10;
      *bound = 1.0e0;
      goto S20;
  S10:
      *bound = 4.0e0;
  S20:
      *status = -1;
      return;
  S30:
      if(*which == 1) goto S70;
  //
  //     P
  //
      if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
      if(!(*p <= 0.0e0)) goto S40;
      *bound = 0.0e0;
      goto S50;
  S40:
      *bound = 1.0e0;
  S50:
      *status = -2;
      return;
  S70:
  S60:
      if(*which == 1) goto S110;
  //
  //     Q
  //
      if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
      if(!(*q <= 0.0e0)) goto S80;
      *bound = 0.0e0;
      goto S90;
  S80:
      *bound = 1.0e0;
  S90:
      *status = -3;
      return;
  S110:
  S100:
      if(*which == 1) goto S150;
  //
  //     P + Q
  //
      pq = *p+*q;
      if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*1.e-16)) goto S140;
      if(!(pq < 0.0e0)) goto S120;
      *bound = 0.0e0;
      goto S130;
  S120:
      *bound = 1.0e0;
  S130:
      *status = 3;
      return;
  S150:
  S140:
      if(*which == 4) goto S170;
  //
  //     SD
  //
      if(!(*sd <= 0.0e0)) goto S160;
      *bound = 0.0e0;
      *status = -6;
      return;
  S170:
  S160:
  //
  //     Calculate ANSWERS
  //
      if(1 == *which) {
  //
  //     Computing P
  //
          z = (*x-*mean)/ *sd;
          cumnor(&z,p,q);
      }
      else if(2 == *which) {
  //
  //     Computing X
  //
          z = dinvnr(p,q);
          *x = *sd*z+*mean;
      }
      else if(3 == *which) {
  //
  //     Computing the MEAN
  //
          z = dinvnr(p,q);
          *mean = *x-*sd*z;
      }
      else if(4 == *which) {
  //
  //     Computing SD
  //
          z = dinvnr(p,q);
          *sd = (*x-*mean)/z;
      }
      return;
  }
  //****************************************************************************
  double dinvnr ( double *p, double *q )

  //****************************************************************************
  // 
  //  Purpose:
  // 
  //    DINVNR computes the inverse of the normal distribution.
  //
  //  Discussion:
  //
  //    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
  //    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
  //
  //    The rational function on page 95 of Kennedy and Gentle is used as a start
  //    value for the Newton method of finding roots.
  //
  //  Reference:
  //
  //    Kennedy and Gentle, 
  //    Statistical Computing,
  //    Marcel Dekker, NY, 1980,
  //    QAA276.4  K46
  //
  //  Parameters:
  //
  //    Input, double *P, *Q, the probability, and the complementary
  //    probability.
  //
  //    Output, double DINVNR, the argument X for which the
  //    Normal CDF has the value P.
  //
  {
  #define maxit 100
  #define eps (1.0e-13)
  #define r2pi 0.3989422804014326e0
  #define nhalf (-0.5e0)
  #define dennor(x) (r2pi*exp(nhalf*(x)*(x)))
  static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
  static int i;
  static unsigned long qporq;

  //
  //     FIND MINIMUM OF P AND Q
  //
      qporq = *p <= *q;
      if(!qporq) goto S10;
      pp = *p;
      goto S20;
  S10:
      pp = *q;
  S20:
  //
  //     INITIALIZATION STEP
  //
      strtx = stvaln(&pp);
      xcur = strtx;
  //
  //     NEWTON INTERATIONS
  //
      for ( i = 1; i <= maxit; i++ )
      {
          cumnor(&xcur,&cum,&ccum);
          dx = (cum-pp)/dennor(xcur);
          xcur -= dx;
          if(fabs(dx/xcur) < eps) goto S40;
      }
      dinvnr = strtx;
  //
  //     IF WE GET HERE, NEWTON HAS FAILED
  //
      if(!qporq) dinvnr = -dinvnr;
      return dinvnr;
  S40:
  //
  //     IF WE GET HERE, NEWTON HAS SUCCEDED
  //
      dinvnr = xcur;
      if(!qporq) dinvnr = -dinvnr;
      return dinvnr;
  #undef maxit
  #undef eps
  #undef r2pi
  #undef nhalf
  #undef dennor
  }
  //****************************************************************************
  double stvaln ( double *p )

  //****************************************************************************
  //
  //  Purpose:
  // 
  //    STVALN provides starting values for the inverse of the normal distribution.
  //
  //  Discussion:
  //
  //    The routine returns X such that 
  //      P = CUMNOR(X),  
  //    that is, 
  //      P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
  //
  //  Reference:
  //
  //    Kennedy and Gentle,
  //    Statistical Computing, 
  //    Marcel Dekker, NY, 1980, page 95,
  //    QAA276.4  K46
  //
  //  Parameters:
  //
  //    Input, double *P, the probability whose normal deviate 
  //    is sought.
  //
  //    Output, double STVALN, the normal deviate whose probability
  //    is P.
  //
  {
  static double xden[5] = {
      0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
      0.38560700634e-2
  };
  static double xnum[5] = {
      -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
      -0.453642210148e-4
  };
  static int K1 = 5;
  static double stvaln,sign,y,z;

      if(!(*p <= 0.5e0)) goto S10;
      sign = -1.0e0;
      z = *p;
      goto S20;
  S10:
      sign = 1.0e0;
      z = 1.0e0-*p;
  S20:
      y = sqrt(-(2.0e0*log(z)));
      stvaln = y+ eval_pol ( xnum, &K1, &y ) / eval_pol ( xden, &K1, &y );
      stvaln = sign*stvaln;
      return stvaln;
  }
  //**********************************************************************
  double eval_pol ( double a[], int *n, double *x )

  //****************************************************************************
  // 
  //  Purpose:
  // 
  //    EVAL_POL evaluates a polynomial at X.
  //
  //  Discussion:
  //
  //    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
  //
  //  Modified:
  //
  //    15 December 1999
  //
  //  Parameters:
  //
  //    Input, double precision A(0:N), coefficients of the polynomial.
  //
  //    Input, int *N, length of A.
  //
  //    Input, double *X, the point at which the polynomial 
  //    is to be evaluated.
  //
  //    Output, double EVAL_POL, the value of the polynomial at X.
  //
  {
    static double devlpl,term;
    static int i;

    term = a[*n-1];
    for ( i = *n-1-1; i >= 0; i-- )
    {
      term = a[i]+term**x;
    }

    devlpl = term;
    return devlpl;
  }
  //****************************************************************************
  
}