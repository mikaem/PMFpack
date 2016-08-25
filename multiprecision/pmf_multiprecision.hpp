#ifndef PMF_MULTIPRECISION_HPP
#define PMF_MULTIPRECISION_HPP

#include <complex>
#include <vector>
#include <algorithm>
#include <limits>
#include <map>
#include <boost/timer/timer.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/function.hpp>
#include <boost/fusion/include/tuple.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "transforms.hpp"

#define NUM_THREADS 2

using std::numeric_limits; 
using boost::timer::auto_cpu_timer;
using boost::timer::cpu_timer;

template <class T>
class DCT
{ // Discrete cosine transform type I. 
  // Direct with high precision, but slow for large N (O(N*N))
public:
  DCT(size_t N)
  : N(N)
  {
    init();
  }
  
  void init()
  {
    static const T pi = boost::math::constants::pi<T>();
    T factor = pi / (N-1);
    for (size_t k=0; k<N; k++){
      for (size_t n=1; n<N-1; n++){
        size_t m = (n*k) % ((N-1)/2);
        auto it = cosmap.find(m);
        if (it == cosmap.end())
          cosmap[m] = cos(factor*m);
      }
    }
  }
  
  std::vector<T> dct_odd(const std::vector<T> &x)
  {
    assert(N % 2 == 1);
    assert(N == x.size());
    
    std::vector<T> xk(N);  
    for (size_t k=0; k<N; k++){
      int sig = (k & 1) ? -1 : 1;
      xk[k] = (x[0]+sig*x[N-1])/2;
      for (size_t n=1; n<N-1; n++){
        T xn = x[n];
        size_t m = ((n*k) % ((N-1)/2));        
        size_t f = ((n*k) % (N-1));      
        size_t q = (n*k) % (2*(N-1));
        
        if (f == 0 and q == 0)
          xk[k] += xn*cosmap[m];
        else if (f == 0 &! q == 0)
          xk[k] -= xn*cosmap[m];
        else if (q > N-1)
        { // lower half
          if (m == 0)
            xk[k] += xn*cosmap[(N-1)/2-m];  
          else if (q < (3*(N-1)/2))
            xk[k] -= xn*cosmap[m];
          else
          xk[k] += xn*cosmap[(N-1)/2-m];
        }
        else
        {
          if (m == 0)
            xk[k] -= xn*cosmap[(N-1)/2-m];  
          else if (q > (N-1)/2)
            xk[k] -= xn*cosmap[(N-1)/2-m];
          else
            xk[k] += xn*cosmap[m];
        }      
      }
    }
    return xk;
  }
  
  T dct_odd_k(size_t k, const std::vector<T> &x)
  {
  //   auto_cpu_timer timer;
    assert(N % 2 == 1);
    assert(N == x.size());
    int sig = (k & 1) ? -1 : 1;
    T xk = (x[0]+sig*x[N-1])/2;
    for (size_t n=1; n<N-1; n++){
      T xn = x[n];
      size_t m = ((n*k) % ((N-1)/2));
      size_t f = ((n*k) % (N-1));      
      size_t q = (n*k) % (2*(N-1));
      
      if (f == 0 and q == 0)
        xk += xn*cosmap[m];
      else if (f == 0 &! q == 0)
        xk -= xn*cosmap[m];
      else if (q > N-1)
      { // lower half
        if (m == 0)
          xk += xn*cosmap[(N-1)/2-m];  
        else if (q < (3*(N-1)/2))
          xk -= xn*cosmap[m];
        else
          xk += xn*cosmap[(N-1)/2-m];
      }
      else
      {
        if (m == 0)
          xk -= xn*cosmap[(N-1)/2-m];  
        else if (q > (N-1)/2)
          xk -= xn*cosmap[(N-1)/2-m];
        else
          xk += xn*cosmap[m];
      }      
    }
    return xk;
  }
    
  size_t N;
  std::map<size_t, T>  cosmap;
};


template <class T>
class DST
{ // Discrete sine transform type I
public:
  DST(size_t N)
  : N(N)
  {
    init();
  }
  
  void init()
  {
    static const T pi = boost::math::constants::pi<T>();
    for (size_t k=0; k<N; k++){
      for (size_t n=0; n<N; n++){
        size_t m = ((n+1)*(k+1)) % (2*(N+1));
        auto it = sinmap.find(m);
        if (it == sinmap.end())
          sinmap[m] = sin(pi*(T(m) / (N+1)));
      }
    }
  }  
  
  std::vector<T> dst_odd(const std::vector<T> &x)
  {
  //   auto_cpu_timer timer;
    assert(N == x.size());
    std::vector<T> xk(N);  
    for (size_t k=0; k<N; k++){
      for (size_t n=0; n<N; n++){
        T xn = x[n];
        size_t m = (((n+1)*(k+1)) % (2*(N+1)));
        xk[k] += xn*sinmap[m];
      }
    }
    return xk;
  }

  T dst_odd_k(size_t k, const std::vector<T> &x)
  {
//     auto_cpu_timer timer;
    assert(N == x.size());
    T xk=0;    
    for (size_t n=0; n<N; n++){
      T xn = x[n];
      size_t m = (((n+1)*(k+1)) % (2*(N+1)));
      xk += xn*sinmap[m];
    }
    return xk;
  }
  
  std::map<size_t, T> sinmap;  
  size_t N;
};

template <class T> 
class ClenshawCurtisBase
{
public:
  
  ClenshawCurtisBase(T epsilon, T fmean, T sigma) 
  : epsilon(epsilon), fmean(fmean), sigma(sigma)
  {
    im = fmean * fmean + sigma;
    read_weights();
  }
  
  ClenshawCurtisBase(T epsilon, T fmean, T sigma, size_t imin, size_t imax) 
  : epsilon(epsilon), fmean(fmean), sigma(sigma), interval_min(imin), interval_max(imax)
  {
    im = fmean * fmean + sigma;
    read_weights();
  }
  
  void read_weights()
  {
    std::ifstream myfile;
    if (std::numeric_limits<T>::digits10 > 32)      
      myfile.open("weights64.dat");
    else
      myfile.open("weights32.dat");
    
    if (myfile.is_open())
    {
      read_weights_from_file(myfile);
      myfile.close();
    }
  }
  
  T integrate_domain(T xmin, T xmax, T tau, T alfa)
  {
    T new_integral = 0;
    T prev_integral = 0;
    std::vector<T> previous_fvals;
    T xmhalf = (xmax - xmin) / 2;
    T xmsum = (xmax + xmin) / 2;    
    
    for (size_t i = interval_min; i < interval_max; i++)
    {
        prev_integral = new_integral;

        size_t Z = pow(2, i) + 1;
        
        fx.resize(Z);
        y.resize(Z);
        
        auto it = weights.find(i);
        if (it == weights.end())
          compute_weights(i);

        std::vector<T> xx = x[i];        
        for (size_t n=0; n<Z; n++)
          y[n] = xx[n]*xmhalf + xmsum;
        
#pragma omp parallel num_threads(NUM_THREADS)
{        
        if (i == this->interval_min)
        {
          #pragma omp for
          for (size_t n=0; n<Z; n++)
            fx[n] = function(y[n], tau, alfa);
        }
        else
        { // Reuse function values from coarser grid
          #pragma omp for
          for (size_t n=1; n<Z; n=n+2)
            fx[n] = function(y[n], tau, alfa);
          #pragma omp for
          for (size_t n=0; n<Z; n=n+2)
            fx[n] = previous_fvals[n/2];
        }
}
        // Store function evaluations for even finer grid
        previous_fvals.resize(Z);
        for (size_t n=0; n<Z; n++)
          previous_fvals[n] = fx[n];
        
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << fx[n] << std::endl;
        
        new_integral = 0;
        
        std::vector<T> ww = weights[i];
        for (size_t n=0; n<Z; n++)
          new_integral += ww[n]*fx[n];   
        
        new_integral = new_integral*xmhalf;

//       std::cout << std::setprecision(numeric_limits<T>::digits10);
//       std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
        if (fabs(new_integral - prev_integral) < epsilon)
          break;        
    }        
    return new_integral;
  }
    
  void set_parameters(T fm, T sm)
  {
    fmean = fm;
    sigma = sm;
    im = fmean * fmean + sigma;    
  }

  void compute_weights(size_t i)
  {
//     auto_cpu_timer timer;  
    size_t Z = pow(2, i) + 1; 
    static const T pi = boost::math::constants::pi<T>();
    size_t M = Z-1;
    std::vector<T> z(Z);
    std::vector<T> w(Z);
    
#pragma omp parallel num_threads(NUM_THREADS)
{
    #pragma omp for
    for (size_t n=0; n<Z; n++)
    {   
      int i = (2*n)-M;       
      z[n] = sin(i*pi/(M*2));
    }
}
    x[i] = z;    
    c.resize(Z);
    
#pragma omp parallel num_threads(NUM_THREADS)
{
    #pragma omp for    
    for (size_t n=0; n<Z; n=n+2)
      c[n] = T(2) / (T(1) - T(n*n));
    
    #pragma omp for    
    for (size_t n=0; n<Z-1; n=n+2)
      c[n+1] = 0;
}    

    DCT<T>* dct = new DCT<T>(Z);
    std::vector<T> ff = dct->dct_odd(c);    
    w[0] = ff[0] / M;
    w[Z-1] = ff[Z-1] / M;
#pragma omp parallel num_threads(NUM_THREADS)
{
    #pragma omp for    
    for (size_t n=1; n<Z-1; n++)
      w[n] = 2*ff[n] / M;
}

//     work.resize(Z*10);    
//     costi(&Z, &work.data()[0], &ifac[0]);
//     cost (&Z, &c.data()[0], &work.data()[0], &ifac[0]);
//     w[0] = c[0] / (2*M);
//     w[Z-1] = c[Z-1] / (2*M);
//     for (size_t n=1; n<Z-1; n++)
//       w[n] = 2*c[n] / (2*M);
        
//     typedef std::complex<T> C;
//     std::vector<C> fxc;
//     fxc.resize(2*M);    
//     for (size_t n=0; n<Z; n++)
//       fxc[n] = c[n];
//         
//     for (size_t n=0; n<Z-2; n++)
//       fxc[Z+n] = c[Z-2-n];
//     
//     fft1(fxc, -1);
//     
//     w[0] = fxc[0].real() / (2*M);
//     w[Z-1] = fxc[Z-1].real() / (2*M);
//     for (size_t n=1; n<Z-1; n++)
//       w[n] = 2*fxc[n].real() / (2*M);
// 
    weights[i] = w;
    
//     std::cout << i << std::endl;
//     for (size_t n=0; n<w.size(); n++)
//     {
//       std::cout << " w "<< n << " " << w[n] << std::endl;
//     }
       
  }
  
  void compute_and_write_weights_to_file(std::ofstream& myfile)
  {
    for (size_t i=interval_min; i<interval_max; i++)
    {
      auto it = weights.find(i);
      if( it == weights.end())
        compute_weights(i);
      
      myfile << std::scientific << std::setprecision(numeric_limits<T>::digits10);
      myfile << weights[i].size() << std::endl;    
      for (size_t n=0; n<weights[i].size(); n++)
        myfile << weights[i][n] << std::endl;
    }
  }
  
  void read_weights_from_file(std::ifstream& myfile)
  {
//     auto_cpu_timer timer;
    static T pi = boost::math::constants::pi<T>();
    std::string line;
    T weight;
    while ( std::getline (myfile,line) )
    {
      int num_weight = boost::lexical_cast<int>(line);
      size_t ii = log2(num_weight-1);
      if (ii > interval_max)
      {
        break;
      }  
      else
      {
        std::vector<T> w(num_weight);
        for (size_t i=0; i<num_weight; i++)
        {
          std::getline (myfile,line);
          weight = boost::lexical_cast<T>(line);
//           std::stringstream ss;          
//           ss << line;
//           ss >> weight;
          
          w[i] = weight;
        }
        weights[ii] = w;
        
        std::vector<T> z(num_weight);    
        size_t M = num_weight-1;
        for (size_t n=0; n<num_weight; n++)
        {   
          int i = n*2-M;
          z[n] = sin(i*pi/(2*M));
        }
        x[ii] = z;
      }
    }   
  }

  virtual T function(T, T, T) = 0;  
  virtual std::vector<T> function_fdf(T, T, T) = 0;
  virtual std::vector<T> integrate_domain_fdf(T xmin, T xmax, T tau, T alfa) = 0;
  virtual T integrate(T) = 0;  
  virtual boost::math::tuple<T, T> integrate_fdf(T tau) = 0;
  T epsilon, sigma, fmean, im, tau;
  size_t interval_min = 2;
  size_t interval_max = 9;
  
protected:
  
  std::map<size_t, std::vector<T> > x;
  std::vector<T> fx, y, c;
  std::map<size_t, std::vector<T> > weights; 
//   int ifac[64];
// std::vector<T> work;

};

template <class T> 
class ClenshawCurtisTransformed: public ClenshawCurtisBase<T>
{
public:
  
  // Declare to avoid the "this" pointer
  using ClenshawCurtisBase<T>::weights;
  using ClenshawCurtisBase<T>::x;
  using ClenshawCurtisBase<T>::fmean;
  using ClenshawCurtisBase<T>::sigma;
  using ClenshawCurtisBase<T>::im;
  using ClenshawCurtisBase<T>::interval_min;
  using ClenshawCurtisBase<T>::interval_max;
  using ClenshawCurtisBase<T>::y;
  using ClenshawCurtisBase<T>::fx;

  ClenshawCurtisTransformed(T epsilon, T fmean, T sigma)
  : ClenshawCurtisBase<T>(epsilon, fmean, sigma)
  {}
  
  ClenshawCurtisTransformed(T epsilon, T fmean, T sigma, size_t imin, size_t imax)
  : ClenshawCurtisBase<T>(epsilon, fmean, sigma, imin, imax)
  {}
    
  T function(T x, T tau, T alfa)
  {
    if (func == 0)
    {
      T d = (1 - x * x);
      T phi = x / d; 
      T X = (1+boost::math::erf((phi - alfa) / (2*sqrt(tau))));
      return X * X * exp(- phi * phi / (2-4*tau)) * (1+x*x) / (d*d);
    }
    else if (func == 1)
    {
      T d = (1 - x);
      T phi = alfa + x / d; 
      T X = (1+boost::math::erf( x / d / (2*sqrt(tau)))) / d;
      return X * X * exp(- phi * phi / (2-4*tau));
    }
    else if (func == 2)
    {
      T d = (1 - x);
      T phi = alfa - d / x; 
      T X = (1+boost::math::erf( -d / x / (2*sqrt(tau)))) / x;
      return X * X * exp(- phi * phi / (2-4*tau));
    }    
  }
  
  T integrate(T tau)
  {
    static const T pi = boost::math::constants::pi<T>();
    T alfa = sqrt(T(2.0)) * boost::math::erfc_inv(2*fmean);
    // Compute constant divisor in integral
    T div = 4*sqrt(pi*(2-4*tau));
    T xmin;
    T xmax;
    if (tau > 0.25)
    {
      xmax = sqrt((2-4*tau)*small2);
      T a = (-1 / xmax + sqrt(1/(xmax*xmax)+4)) / 2;
      T b = (-1 / xmax - sqrt(1/(xmax*xmax)+4)) / 2;
      xmax = a*a < 1 ? a : b;             
      
      xmin = (alfa - 2*sqrt(tau)*small);  
      a = (-1 / xmin + sqrt(1/(xmin*xmin)+4)) / 2;
      b = (-1 / xmin - sqrt(1/(xmin*xmin)+4)) / 2;
      xmin = a*a < 1 ? a : b;             
//         std::cout << xmin << " " << xmax <<  " " << a << " " << b << std::endl;      
      func = 0;
      integral = this->integrate_domain(xmin, xmax, tau, alfa); 
    }
    else
    {    
      xmin = T(1) / (1+2*sqrt(tau)*small);
      func = 2;
      integral = this->integrate_domain(xmin, T(1), tau, alfa); 
      
      xmax = sqrt((2-4*tau)*small2);
      xmax = (xmax-alfa) / (1+xmax-alfa);      
//         std::cout << xmin << " " << xmax << std::endl;
      func = 1;
      integral  += this->integrate_domain(T(0), xmax, tau, alfa); 
    }
    
//      std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
//      std::cout <<  tau << " " << T(1.0) - integral / im / div << std::endl;

    return T(1) - integral / im / div;    
  }
  
  std::vector<T> function_fdf(T x, T tau, T alfa)
  {
    static const T root_pi = boost::math::constants::root_pi<T>();
    if (func == 0)
    {
      T d = (1 - x * x);
      T phi = x / d; 
      T tmp = (phi - alfa) / (2*sqrt(tau));
      T X = (1+boost::math::erf(tmp));
      T df = X * (alfa - phi / (1-2*tau)) * exp(-tmp*tmp-phi*phi/(2-4*tau));
      std::vector<T> result = {X * X * exp(- phi * phi / (2-4*tau)) * (1+x*x) / (d*d),
        df * (1+x*x) / (d*d)};
      return result;
    }
    else if (func == 1)
    {
      T d = (1 - x);
      T phi = alfa + x / d; 
      T tmp = x / d / (2*sqrt(tau));
      T X = (1+boost::math::erf(tmp));
      T df = X * (alfa - phi / (1-2*tau)) * exp(-tmp*tmp-phi*phi/(2-4*tau));
      std::vector<T> result = {X * X * exp(- phi * phi / (2-4*tau)) / (d*d),
        df / (d*d)};
      return result;
    }
    else if (func == 2)
    {
      T d = (1 - x);
      T phi = alfa - d / x; 
      T tmp = - d / x / (2*sqrt(tau));
      T X = (1+boost::math::erf(tmp));
      T df = X * (alfa - phi / (1-2*tau)) * exp(-tmp*tmp-phi*phi/(2-4*tau));
      std::vector<T> result = {X * X * exp(- phi * phi / (2-4*tau)) / (x*x),
        df / (x*x)};
      return result;
    }    
  }

  boost::math::tuple<T, T> integrate_fdf(T tau)
  {
    static const T pi = boost::math::constants::pi<T>();
    static const T root_pi = boost::math::constants::root_pi<T>();
    T alfa = sqrt(T(2)) * boost::math::erfc_inv(2*fmean);
    // Compute constant divisor in integral
    T div = 4*sqrt(pi*(2-4*tau));
    T xmin;
    T xmax;
    std::vector<T> integrals;
    
    if (tau > 0.25)
    {
      xmax = sqrt(2*(1-2*tau)*small2);
      T a = (-1 / xmax + sqrt(1/(xmax*xmax)+4)) / 2;
      T b = (-1 / xmax - sqrt(1/(xmax*xmax)+4)) / 2;
      xmax = a*a < 1 ? a : b;             
      
      xmin = (alfa - 2*sqrt(tau)*small);  
      a = (-1 / xmin + sqrt(1/(xmin*xmin)+4)) / 2;
      b = (-1 / xmin - sqrt(1/(xmin*xmin)+4)) / 2;
      xmin = a*a < 1 ? a : b;             
      
      func = 0;
      integrals = this->integrate_domain_fdf(xmin, xmax, tau, alfa); 
    }
    else
    {    
      xmin = T(1) / (1+2*sqrt(tau)*small);
      func = 2;
      integrals = this->integrate_domain_fdf(xmin, T(1), tau, alfa); 
            
      xmax = sqrt((2-4*tau)*small2)-alfa;
      xmax = xmax / (1+xmax);      
      
      func = 1;
      std::vector<T> fdintegrals  = this->integrate_domain_fdf(T(0), xmax, tau, alfa); 
      integrals[0] += fdintegrals[0];
      integrals[1] += fdintegrals[1];
    }
    
    // Scale the integrals appropriately
    integrals[1] = integrals[1] / tau / sqrt(tau) / root_pi;
    integrals[1] = -integrals[1] / im / div;
    
    integrals[0] = T(1) - integrals[0] / im / div;

    return boost::math::make_tuple(integrals[0], integrals[1]);    
  }
  
  std::vector<T> integrate_domain_fdf(T xmin, T xmax, T tau, T alfa)
  {
    T new_integral = 0;
    T fdf_integral = 0;
    T prev_integral = 0;
    std::vector<T> previous_fvals;
    std::vector<T> previous_fdfvals;
    T xmhalf = (xmax - xmin) / 2;
    T xmsum = (xmax + xmin) / 2;  
    
    for (size_t i = interval_min; i < interval_max; i++)
    {
      prev_integral = new_integral;

      size_t Z = pow(2, i) + 1;
      
      fdf.resize(Z);
      fx.resize(Z);
      y.resize(Z);
      
      auto it = weights.find(i);
      if (it == weights.end())
        this->compute_weights(i);

      auto xx = x[i];
#pragma omp parallel num_threads(NUM_THREADS)
{    
      #pragma omp for
      for (size_t n=0; n<Z; n++)
        y[n] = xx[n]*xmhalf + xmsum;
      
      if (i == interval_min)
      {
        #pragma omp for
        for (size_t n=0; n<Z; n++)
        {
          auto vals = function_fdf(y[n], tau, alfa);
          fx[n] = vals[0];
          fdf[n] = vals[1];
        }
      }
      else
      { // Reuse function values from coarser grid
        #pragma omp for
        for (size_t n=1; n<Z; n=n+2)
        {
          auto vals = function_fdf(y[n], tau, alfa);
          fx[n] = vals[0];
          fdf[n] = vals[1];
        }
        #pragma omp for
        for (size_t n=0; n<Z; n=n+2)
        {
          fx[n] = previous_fvals[n/2];
          fdf[n] = previous_fdfvals[n/2];
        }
      }
}      
      // Store function evaluations for even finer grid
      previous_fvals.resize(Z);
      previous_fdfvals.resize(Z);
      for (size_t n=0; n<Z; n++)
      {
        previous_fvals[n] = fx[n];
        previous_fdfvals[n] = fdf[n];
      }
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << fx[n] << " " << fdf[n] << std::endl;
      
      new_integral = T(0);
      std::vector<T> ww = weights[i];
            
      for (size_t n=0; n<Z; n++)
        new_integral += ww[n]*fx[n];   
      
      new_integral = new_integral*xmhalf;

      fdf_integral = T(0);
      for (size_t n=0; n<Z; n++)
        fdf_integral += ww[n]*fdf[n];   
      
      fdf_integral = fdf_integral*xmhalf;
      
//        std::cout << std::setprecision(numeric_limits<T>::digits10);
//        std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
      if (fabs(new_integral - prev_integral) < this->epsilon)
        break;        
    }
    std::vector<T> result = {new_integral, fdf_integral};
    return result;
  }
  
  std::vector<T> fdf;
  T integral;
  T small = boost::math::erfc_inv(this->epsilon/100);
//   T small2 = -log(this->epsilon/10000000);
  T small2 = 10*std::numeric_limits<T>::digits10;
  size_t func;
};

template <class T> 
class ClenshawCurtis: public ClenshawCurtisBase<T>
{
public:
  
  using ClenshawCurtisBase<T>::fx;
  using ClenshawCurtisBase<T>::y;
  
  ClenshawCurtis(T epsilon, T fmean, T sigma)
  : ClenshawCurtisBase<T>(epsilon, fmean, sigma)
  {}
  
  ClenshawCurtis(T epsilon, T fmean, T sigma, size_t imin, size_t imax)
  : ClenshawCurtisBase<T>(epsilon, fmean, sigma, imin, imax)
  {}
    
  T function(T x, T tau, T alfa)
  {
    T X = (1+boost::math::erf((x - alfa) / (2*sqrt(tau))));
    return X * X * exp(- x * x / (2-4*tau));
  }
  
  std::vector<T> function_fdf(T x, T tau, T alfa)
  {
     T tmp = (x - alfa) / (2*sqrt(tau));
     T X = (1+boost::math::erf(tmp));
     T f = X * X * exp(- x * x / (2-4*tau));
     T df = X * (alfa - x / (1-2*tau)) * exp(-tmp*tmp-x*x/(2-4*tau));
     std::vector<T> result = {f, df};
     return result;
  }

  T integrate(T tau)
  {
    static const T pi = boost::math::constants::pi<T>();
    T alfa = sqrt(T(2)) * boost::math::erfc_inv(2*this->fmean);
    T integral = 0;
    T integral_dom = 1;
    T two_root_tau = 2*sqrt(tau);
    T phi_X_unity = (alfa + two_root_tau * this->inv_eps) ; // above this X in integrand is > 1-epsilon and can be set to unity, which makes the integral analytical
    T phi_X_99 = (alfa + two_root_tau * this->inv_99);
    T phi_X_01 = (alfa + two_root_tau * this->inv_01);
    T phi_X_5 = alfa;
    
    // Compute constant divisor in integral
    T div = 4*sqrt(pi*(2-4*tau));
    
    // For low variance the integral is smoother, no sharp gradient around alfa
    if (tau > 0.25)
    {
      phi_X_5 = 0;
      phi_X_99 = sqrt((4*tau-2)*log(div / 400));
      phi_X_01 = -phi_X_99;
    }
    T xmin = phi_X_5;
    T xmax = 0;
    size_t i = 1;
    
    T dphi = (phi_X_99-phi_X_01) / 2;
    dphi = dphi < T(2) / 10 ? dphi : T(2) / 10;
    right = true;
    
    while (integral_dom / integral > this->epsilon*this->sigma)
    { 
       xmax = xmin + T(i) * dphi;
       
       integral_dom = this->integrate_domain(xmin, xmax, tau, alfa);
       integral_dom /= div;
       integral += integral_dom;
       xmin = xmax;
       i++;
              
       // Check if xmin is large such that X in the integrand is unity
       if (xmin > phi_X_unity)
       {
         integral += (1 - cdf(this->s, xmin/sqrt(1-2*tau)));
         break;
       }        
    }
    
    // Now move to the left of phi_X_5
    this->f_lower = 0;
    this->f_upper = 0;
    xmax = phi_X_5;
    integral_dom = 1;
    i = 1;
    right = false;
    while (integral_dom > this->epsilon*this->sigma)
    {
       xmin = xmax - T(i) * dphi;

       integral_dom = this->integrate_domain(xmin, xmax, tau, alfa);
       integral_dom /= div;
       integral += integral_dom;
       xmax = xmin;
       i++;
    }
    
//     std::cout <<  tau << " " << T(1.0) - integral / this->im << std::endl;
    
    return T(1) - integral / this->im;    
  }
  
  T integrate_domain(T xmin, T xmax, T tau, T alfa)
  {
    T new_integral = 0;
    T prev_integral = 0;
    T f_shared = 0;
    std::vector<T> previous_fvals;
    T xmhalf = (xmax - xmin) / 2;
    T xmsum = (xmax + xmin) / 2;    
    
    for (size_t i = this->interval_min; i < this->interval_max; i++)
    {
        prev_integral = new_integral;

        size_t Z = pow(2, i) + 1;
        
        fx.resize(Z);
        y.resize(Z);
        
        auto it = this->weights.find(i);
        if (it == this->weights.end())
          this->compute_weights(i);

        std::vector<T> xx = this->x[i];
#pragma omp parallel for num_threads(NUM_THREADS)        
        for (size_t n=0; n<Z; n++)
          y[n] = xx[n]*xmhalf + xmsum;
        
        if (i == this->interval_min)
        {
          if (right and f_upper > this->epsilon)
          {
            fx[0] = f_upper;
#pragma omp parallel for num_threads(NUM_THREADS)        
            for (size_t n=1; n<Z; n++)
              fx[n] = function(y[n], tau, alfa);
          }
          else if (not right and f_lower > this->epsilon)
          {
            fx[Z-1] = f_lower;
#pragma omp parallel for num_threads(NUM_THREADS)        
            for (size_t n=0; n<Z-1; n++)
              fx[n] = function(y[n], tau, alfa);
          }
          else
          {           
#pragma omp parallel for num_threads(NUM_THREADS)        
            for (size_t n=0; n<Z; n++)
              fx[n] = function(y[n], tau, alfa);
          }
        }
        else
        { // Reuse function values from coarser grid
#pragma omp parallel for num_threads(NUM_THREADS)        
          for (size_t n=1; n<Z; n=n+2)
            fx[n] = function(y[n], tau, alfa);
#pragma omp parallel for num_threads(NUM_THREADS)        
          for (size_t n=0; n<Z; n=n+2)
            fx[n] = previous_fvals[n/2];
        }
        // Store function evaluations for even finer grid
        previous_fvals.resize(Z);
#pragma omp parallel for num_threads(NUM_THREADS)        
        for (size_t n=0; n<Z; n++)
          previous_fvals[n] = fx[n];
        
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << fx[n] << std::endl;
        
        // Store edge values for next interval
        f_upper = fx[Z-1];
        f_lower = fx[0];
         
        new_integral = 0;
        
        std::vector<T> ww = this->weights[i];
        for (size_t n=0; n<Z; n++)
          new_integral = new_integral + ww[n]*fx[n];   
        new_integral = new_integral*xmhalf;
        
//          std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
        if (fabs(new_integral - prev_integral) < this->epsilon)
          break;        
    }        
    return new_integral;
  }

  std::vector<T> integrate_domain_fdf(T xmin, T xmax, T tau, T alfa)
  {
    T new_integral = 0;
    T prev_integral = 0;
    T fdf_integral = 0;
    std::vector<T> previous_fvals;
    std::vector<T> previous_fdfvals;
    T xmhalf = (xmax - xmin) / 2;
    T xmsum = (xmax + xmin) / 2;    
    std::vector<T> vals;
    
    for (size_t i = this->interval_min; i < this->interval_max; i++)
    {
        prev_integral = new_integral;

        size_t Z = pow(2, i) + 1;
        
        fx.resize(Z);
        this->y.resize(Z);
        fdf.resize(Z);
        
        auto it = this->weights.find(i);
        if (it == this->weights.end())
          this->compute_weights(i);

        std::vector<T> xx = this->x[i];
        for (size_t n=0; n<Z; n++)
          this->y[n] = xx[n]*xmhalf + xmsum;
        
        if (i == this->interval_min)
        {
          if (right and f_upper > this->epsilon)
          {
            fx[0] = f_upper;
            fdf[0] = fdf_upper;
            for (size_t n=1; n<Z; n++)
            {
              vals = function_fdf(this->y[n], tau, alfa);
              fx[n] = vals[0];
              fdf[n] = vals[1];
            }
          }
          else if (not right and f_lower > this->epsilon)
          {
            fx[Z-1] = f_lower;
            fdf[Z-1] = fdf_lower;
            for (size_t n=0; n<Z-1; n++)
            {
               vals = function_fdf(this->y[n], tau, alfa);
               fx[n] = vals[0];
               fdf[n] = vals[1];              
            }
          }
          else
          {           
            for (size_t n=0; n<Z; n++)
            {
              vals = function_fdf(this->y[n], tau, alfa);
              fx[n] = vals[0];
              fdf[n] = vals[1];
            }
          }
        }
        else
        { // Reuse function values from coarser grid
          for (size_t n=1; n<Z; n=n+2)
          {
            vals = function_fdf(this->y[n], tau, alfa);
            fx[n] = vals[0];
            fdf[n] = vals[1];
          }
          for (size_t n=0; n<Z; n=n+2)
          {
            fx[n] = previous_fvals[n/2];
            fdf[n] = previous_fdfvals[n/2];
          }          
        }
        // Store function evaluations for even finer grid
        previous_fvals.resize(Z);
        previous_fdfvals.resize(Z);
        for (size_t n=0; n<Z; n++)
        {
          previous_fvals[n] = fx[n];
          previous_fdfvals[n] = fdf[n];
        }
        
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << fx[n] << std::endl;
        
        // Store edge values for next interval
        f_upper = fx[Z-1];
        f_lower = fx[0];
        fdf_upper = fdf[Z-1];
        fdf_lower = fdf[0];
         
        new_integral = 0;
        
        std::vector<T> ww = this->weights[i];
        for (size_t n=0; n<Z; n++)
          new_integral = new_integral + ww[n]*fx[n];   
        new_integral = new_integral*xmhalf;
        
        fdf_integral = 0;
        for (size_t n=0; n<Z; n++)
          fdf_integral = fdf_integral + ww[n]*fdf[n];   
        fdf_integral = fdf_integral*xmhalf;
        
//          std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
        if (fabs(new_integral - prev_integral) < this->epsilon)
          break;        
    }        
    std::vector<T> result = {new_integral, fdf_integral};
    return result;

  }
  
  boost::math::tuple<T, T> integrate_fdf(T tau)
  {
    static const T pi = boost::math::constants::pi<T>();
    static const T root_pi = boost::math::constants::root_pi<T>();
    T alfa = sqrt(T(2)) * boost::math::erfc_inv(2*this->fmean);
    T integral = 0;
    T integral_dom = 1;
    T two_root_tau = 2*sqrt(tau);
    T phi_X_unity = (alfa + two_root_tau * this->inv_eps) ; // above this X in integrand is > 1-epsilon and can be set to unity, which makes the integral analytical
    T phi_X_99 = (alfa + two_root_tau * this->inv_99);
    T phi_X_01 = (alfa + two_root_tau * this->inv_01);
    T phi_X_5 = alfa;
    std::vector<T> integrals = {0, 0};
    std::vector<T> integrals_dom = {1, 0};
    
    // Compute constant divisor in integral
    T div = 4*sqrt(pi*(2-4*tau));
    
    // For low variance the integral is smoother, no sharp gradient around alfa
    if (tau > 0.25)
    {
      phi_X_5 = 0;
      phi_X_99 = sqrt((4*tau-2)*log(div / 400));
      phi_X_01 = -phi_X_99;
    }
    T xmin = phi_X_5;
    T xmax = 0;
    size_t i = 1;
    
    T dphi = (phi_X_99-phi_X_01) / 2;
    dphi = dphi < T(2) / 10 ? dphi : T(2) / 10;
    right = true;
    
    while (integrals_dom[0] / integrals[0] > this->epsilon*this->sigma)
    { 
       xmax = xmin + T(i) * dphi;
       
       integrals_dom = this->integrate_domain_fdf(xmin, xmax, tau, alfa);
       integrals_dom[0] /= div;
       integrals_dom[1] /= div;
       integrals[0] += integrals_dom[0];
       integrals[1] += integrals_dom[1];
       xmin = xmax;
       i++;
              
       // Check if xmin is large such that X in the integrand is unity
       if (xmin > phi_X_unity)
       {
         T small2 = 6*std::numeric_limits<T>::digits;
         xmax = sqrt((2-4*tau)*small2);
         integrals_dom = this->integrate_domain_fdf(xmin, xmax, tau, alfa);
         integrals[0] += (1 - cdf(this->s, xmin/sqrt(1-2*tau))) ;
         integrals[1] += integrals_dom[1] / div;
         break;
       }        
    }
    
    
    // Now move to the left of phi_X_5
    this->f_lower = 0;
    this->f_upper = 0;
    fdf_lower = 0;
    fdf_upper= 0;
    xmax = phi_X_5;
    integral_dom = 1;
    integrals_dom[0] = 1;
    i = 1;
    right = false;
    while (integrals_dom[0] > this->epsilon*this->sigma)
    {
       xmin = xmax - T(i) * dphi;

       integrals_dom = this->integrate_domain_fdf(xmin, xmax, tau, alfa);
       integrals_dom[0] /= div;
       integrals_dom[1] /= div;
       integrals[0] += integrals_dom[0];
       integrals[1] += integrals_dom[1];
       xmax = xmin;
       i++;
    }
    // Scale the integrals appropriately
    integrals[1] = integrals[1] / tau / sqrt(tau) / root_pi;
    integrals[1] = -integrals[1] / this->im;
    
    integrals[0] = T(1) - integrals[0] / this->im;
    
//     std::cout <<  tau << " " << T(1.0) - integral / this->im << std::endl;
    
    return boost::math::make_tuple(integrals[0], integrals[1]);        
  }

    
  boost::math::normal_distribution<T> s;
  T inv_eps = -quantile(s, numeric_limits<T>::epsilon()/1000);
  T inv_99  = -quantile(s, T(1) / 100);
  T inv_01  = -quantile(complement(s, T(1) / 100));
  bool right;
  T f_upper, f_lower;
  T fdf_upper, fdf_lower;
  std::vector<T> fdf;
};

template <class T>
struct BoostFunction
{
  BoostFunction(ClenshawCurtisBase<T> *cc) 
  : cc(cc) {};
  T operator ()(T const& z) 
    {
       T tmp = cc->integrate(z); 
//        std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
//        std::cout << tmp  << std::endl;      
       return tmp;
    };
private:
  ClenshawCurtisBase<T> *cc;
};

template <class T>
struct Boost_df_Function
{
  Boost_df_Function(ClenshawCurtisBase<T> *cc) 
  : cc(cc) {};  
  boost::math::tuple<T, T> operator ()(T const& z) 
    {
      boost::math::tuple<T, T> tmp = cc->integrate_fdf(z);
//         std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
//         std::cout << " " << std::get<0>(tmp) << " " << std::get<1>(tmp) << std::endl;      
      return tmp;
    };
private:
  ClenshawCurtisBase<T> *cc;
};

template <class T>
class tolerance {
public:
  tolerance(T eps) :
    _eps(eps) {
  }
  bool operator()(T a, T b) {
    return (fabs(b - a) <= _eps);
  }
private:
  T _eps;
};

template <class T>
class ComputeTau
{
public:
  ComputeTau(T fm, T sm)
  : fmean(fm), sigma(sm)
  {
    epsilon = std::numeric_limits<T>::epsilon();
    digits = std::numeric_limits<T>::digits10;
    tol = new tolerance<T>(epsilon);
    init();
  }

  ComputeTau(T fm, T sm, size_t imin, size_t imax)
  :  fmean(fm), sigma(sm), min_integration_interval(imin), 
     max_integration_interval(imax)
  {
    epsilon = std::numeric_limits<T>::epsilon();
    digits = std::numeric_limits<T>::digits10;
    tol = new tolerance<T>(epsilon);
    init();
  }
  
  ComputeTau(T fm, T sm, size_t imin, size_t imax, T eps, T digs)
  :  fmean(fm), sigma(sm), min_integration_interval(imin), 
     max_integration_interval(imax), epsilon(eps), digits(digs)
  {
    tol = new tolerance<T>(epsilon);
    init();
  }
  
  ~ComputeTau()
  {
    delete CC;
    delete f;
    delete fdf;
    delete tol;
  }
  
  bool set_parameters(T fm, T sm)
  {
    fmean = fm;
    sigma = sm;
    T seg = sigma / fmean / (1-fmean);
    if (seg > (1-std::numeric_limits<T>::epsilon()) or 
      seg < std::numeric_limits<T>::epsilon())
      return false;
    else
    {
      CC->set_parameters(fm, sm);    
      return true;
    }    
  }
  
  void init()
  { 
    CC = new ClenshawCurtisTransformed<T>(epsilon, fmean, sigma, 
                               min_integration_interval, 
                               max_integration_interval);
    f = new BoostFunction<T>(CC);
    fdf = new Boost_df_Function<T>(CC);    
  }

  T compute_newton(T guess, boost::uintmax_t iters=20)
  {
    T tmax = guess*1.05;
    T tmin = guess*0.95;
    T res = boost::math::tools::newton_raphson_iterate(*fdf, guess, tmin, tmax, digits, iters);
//      std::cout << CC->integrate(res) << std::endl;
    return res;
  }
  
  T compute_bracket(T lower=T(1)/100, T upper=T(5)/10-T(1)/100, T beta=T(1)/10)
  { // Find tau using robust bracketing algorithm
    // Use lower and upper as initial guess for bracket. 0 < beta < 1 decides how fast 
    // to reduce or increase lower/upper in effort to bracket root initially.  
    boost::uintmax_t it = 100;  

    T fa = f->operator()(lower);
    T fb = f->operator()(upper);
    int count = 0;
    while (fa > 0 and count < 12)
    { // With beta=0.01 and lower=0.01 -> 0.01, 0.001, 0.0001, ...
      lower = lower * beta;
      fa = f->operator()(lower);
//       std::cout <<fa << " " << lower << std::endl;
      count++;
    } 
        
    count = 0;
    while (fb < 0 and count < 12)
    {  // With beta=0.01 and upper=0.49 -> 0.49, 0.499, 0.4999, ...
      upper = (1-beta) / 2 + upper * beta;
      fb = f->operator()(upper);
//       std::cout <<fb << " " << upper << std::endl;
      count++;
    }
//     std::cout << fa << " " << fb << " " << lower << " " << upper << std::endl;
    
    if (fa*fb < 0)
    {
      std::pair<T, T> found = boost::math::tools::toms748_solve(*f, lower, upper, fa, fb, *tol, it);
      if (it == 100) std::cout << " Root not found Boost toms748 Fsolver" << std::endl;
      tau = (found.first + found.second) / 2;
//     std::cout << tau << std::endl;
      return tau;
    }
    else
    {
      std::cout << " Root not found, could not bracket root. " << std::endl;
      std::cout << fmean << " " << sigma <<  " " << sigma / fmean / (1-fmean) << std::endl;
      std::cout << fa << " " << fb << std::endl;
//       std::cout << lower << " " << upper << std::endl;
      return -1;
    }
  }
    
  ClenshawCurtisTransformed<T>* CC;
  BoostFunction<T>* f;
  Boost_df_Function<T>* fdf;
  T epsilon;
  tolerance<T>* tol;
  T fmean, sigma, tau;
  int digits;
  int mode;
  size_t min_integration_interval = 2;
  size_t max_integration_interval = 9;
};

template <class T, class M>
struct ComputeTauBracketInit
{
  ComputeTauBracketInit(T fmean, T sigma)
  {
    init(fmean, sigma);
  };  

  ComputeTauBracketInit(T fmean, T sigma, size_t imin, size_t imax)
  : min_integration_interval(imin), max_integration_interval(imax)
  {
    init(fmean, sigma);
  };  
  
  void init(T fmean, T sigma)
  {
    ctau = new ComputeTau<T>(fmean, sigma, min_integration_interval, max_integration_interval);
    ftau = new ComputeTau<M>(M(fmean), M(sigma), 2, 7);
  }
  
  ~ComputeTauBracketInit()
  {
    delete ctau;
    delete ftau;
  }
  
  T compute(T fm, T sm) 
  {
    bool ok = ctau->set_parameters(fm, sm);
    ok = ftau->set_parameters(M(fm), M(sm));
    if (ok)
    {
      M tau0 = ftau->compute_bracket();      
      return ctau->compute_newton(T(tau0));
    }
    else
    {
      return -1;
    }
  };
  
  T compute() 
  {
    M tau0 = ftau->compute_bracket();      
    return ctau->compute_newton(T(tau0));
  };
  
private:
  ComputeTau<T>* ctau;
  ComputeTau<M>* ftau;
  size_t mode;
  M tau0;
  size_t min_integration_interval = 2;
  size_t max_integration_interval = 9;
};

template <class T>
class DerivativesHP
{
  // Compute a quick initial guess using bracketing float solver
  // Compute refinement to 1e-16 using newton with double precision
  // Polish root using full precision

public:
  
  DerivativesHP(T fmean, T sigma, T dx, size_t order, bool central, 
                size_t imin=2, size_t imax=9, size_t max_iters_hp=2)
  : fmean(fmean), sigma(sigma), dx(dx), order(order), central(central),  
    min_integration_interval(imin), max_integration_interval(imax),
    max_iters_hp(max_iters_hp)
  {
     dst = new DST<T>(order-2);
     dct = new DCT<T>(order);
  };
    
  ~DerivativesHP()
  {
    delete ctau;
    delete ftau;
    delete dtau;
    delete dst;
    delete dct;
  }
  
  virtual void init()
  {
    ctau = new ComputeTau<T>(fmean, sigma, 
                              min_integration_interval, 
                              max_integration_interval,
                              std::numeric_limits<T>::epsilon(), 
                              2*std::numeric_limits<T>::digits10);
    
    dtau = new ComputeTau<double>(double(fmean), double(sigma), 2, 9, 
                               2*std::numeric_limits<double>::epsilon(), 
                               2*std::numeric_limits<double>::digits10);
    
    ftau = new ComputeTau<float>(float(fmean), float(sigma), 2, 8,
                               5*std::numeric_limits<float>::epsilon(), 
                               std::numeric_limits<float>::digits10);
    
  }
  
  virtual bool set_parameters(T fm, T sm)
  {
    if (!ctau)
      init();

    fmean = fm;
    sigma = sm;
    T seg = sm / fm / (1-fm);
    if (seg > (1-std::numeric_limits<T>::epsilon()) or
        seg < std::numeric_limits<T>::epsilon())
    {
      std::cout << "fmean " << fm << " sigma " << sm << " Seg " << seg << std::endl;
      return false;
    }
    else
    {
      bool ok;
      ok = ctau->set_parameters(fmean, sigma);
      ok = ftau->set_parameters(float(fmean), float(sigma));
      ok = dtau->set_parameters(double(fmean), double(sigma));
      return ok;
    }
  }
  
  virtual T evaluate(T const z) 
  {
    T fm0 = fmean;
    T sigma0 = sigma;
    T im0 = fmean * fmean + sigma;
    if (mode == 0)
    { // Derivative with respect to fmean (sigma constant)
      fmean = z;
    }
    else if (mode == 1)
    { // Derivative with respect to fmean (im constant)
      fmean = z;
      sigma = im0 - fmean * fmean;
    }
    else if(mode == 2)
    { // Derivative with respect to sigma (fmean constant)
      sigma = z;
    }
    else if(mode == 3)
    {  // Derivative with respect to im (fmean constant)
      sigma = z - fmean * fmean;
    }
//     std::cout << fmean << " " << sigma << " " << sigma/fmean/(1-fmean) << std::endl;

    bool ok = set_parameters(fmean, sigma);
    if (not ok)
    {
      return -1;
    }
    else
    {    
      T seg = sigma / fmean / (1-fmean);
      T result;
      
      if (seg > 0.999)
      { // Very high segregation - double gets into trouble due to roundoff. Use robust bracket with HP.
        result = ctau->compute_bracket(T(1)/10000000, T(1)/100);  
      }      
      else if (seg > 0.95)
      { // float cannot find tau anymore due to roundoff. Use robust bracketing with double.
        tau1 = dtau->compute_bracket(0.0001, 0.1);
        result = ctau->compute_newton(T(tau1), max_iters_hp);
      }
      else 
      { // Find first estimate using float (accuracy ~1e-7). Improve accuracy
        // to ~1e-15 using double and then finally a few iterations with the
        // (expensive) high precision solver.
        tau0 = ftau->compute_bracket();
        tau1 = dtau->compute_newton(double(tau0), 3);       
        result = ctau->compute_newton(T(tau1), max_iters_hp);
      }
            
      ctau->set_parameters(fm0, sigma0);
      fmean = fm0;
      sigma = sigma0;
      
      return result;
    }
  };
  
  virtual T compute_tau(T fmean, T sigma)
  {
    set_parameters(fmean, sigma);
    double tau1 = dtau->compute_newton(double(tau0), 3);     
    return ctau->compute_newton(T(tau1), max_iters_hp);
  }
  
  virtual T compute_tau()
  {
    if (!ctau)
      init();

    double tau1 = dtau->compute_newton(double(tau0), 3);     
    return ctau->compute_newton(T(tau1), max_iters_hp);
  }
    
  virtual std::vector<T> derivatives()
  {
//     auto_cpu_timer timer;  

    const T half = boost::math::constants::half<T>();
    const T pi = boost::math::constants::pi<T>();
        
    x = fmean;
    mode = central ? 0 : 1;

    int N = order - 1;
    std::vector<T> z(order);
    std::vector<T> f(order);
    f[N/2] = evaluate(x);
//     T h = dx * f[N/2];   
    T h = dx;   
    T bma = h * half;
    T con = 1 / bma;
    
    for (size_t n=0; n<order; n++)
    {
      if (n != N/2)
      { 
        z[n] = cos(n*pi/N);
        f[n] = evaluate(z[n]*bma + x);
        if (f[n] < 0)
        {
          std::cerr << "Intensity of segregation out of range. Use smaller domain for computing derivatives!!" << std::endl;
        }
      }

//         std::cout << f[n] << " " << z[n]*bma + x << std::endl;
    }

//     for (size_t i=0; i< f.size(); i++)
//       std::cout << f[i] << std::endl;
    
    std::vector<T> dcf = dct->dct_odd(f);
    for (size_t i=0; i<dcf.size() ; i++)
      dcf[i] = dcf[i] * 2;
    
//     for (size_t i=0; i< dcf.size(); i++)
//       std::cout << dcf[i] << std::endl;
    
    std::vector<T> dfs;
    std::vector<T> d2fs;
    d2fs.push_back(0);
    for (int n=1; n<N; n++)
    {      
      dfs.push_back(-dcf[n]*n);
      d2fs.push_back(dcf[n]*n*n);
    }
    d2fs.push_back(dcf[N]*N*N);
    
    std::vector<T> dffs = dst->dst_odd(dfs);
    std::vector<T> d2ffs = dct->dct_odd(d2fs);
    for (size_t i=0; i< dffs.size(); i++)
      dffs[i] = dffs[i] * 2;
    for (size_t i=0; i< d2ffs.size(); i++)
      d2ffs[i] = d2ffs[i] * 2;
            
    std::vector<T> result {f[order/2],
                          -dffs[(order-2)/2] / (2*N) * con, 
                          -d2ffs[order/2] / (2*N) * con * con,
                          0, 0};    
                              
    mode = central ? 2 : 3;
    x = central ? sigma : fmean * fmean + sigma;
    
    for (size_t n=0; n<order; n++)
    {
      if (n != N/2)
        f[n] = evaluate(z[n]*bma + x);
    }

    dcf = dct->dct_odd(f);
    
    for (size_t i=0; i<dcf.size() ; i++)
      dcf[i] = dcf[i] * 2;

    for (int n=1; n<N; n++)
    {      
      dfs[n-1] = -dcf[n]*n;
      d2fs[n]  = dcf[n]*n*n;
    }
    d2fs[N] = dcf[N]*N*N;
    
    T dffs1 = dst->dst_odd_k((order-2)/2, dfs);
    T d2ffs1 = dct->dct_odd_k(order/2, d2fs);
        
    result[3] = -dffs1 / N * con;
    result[4] = -d2ffs1 / N * con * con;
    
    return result;
// //     auto_cpu_timer timer;  
//     typedef std::complex<T> C;  
//     std::vector<C> df;
//     std::vector<C> df2;
//     std::vector<C> dd;
//     
//     if (tau0 < 0)
//       init();
// 
//     const T half = boost::math::constants::half<T>();
//     const T pi = boost::math::constants::pi<T>();
//     T bma = dx * half;
//     T con = 1 / bma;
//     
//     x = fmean;
//     mode = central ? 0 : 1;
// 
//     int N = order - 1;
//     std::vector<T> z(order);
//     std::vector<T> f(order);
//     for (size_t n=0; n<order; n++)
//     {
//       z[n] = cos(n*pi/N);
//       f[n] = evaluate(z[n]*bma + x);
//     }
// 
//     dd.resize(2*N);    
//     for (size_t n=0; n<order; n++)
//       dd[n] = f[n];
//         
//     for (size_t n=0; n<order-2; n++)
//       dd[order+n] = f[order-2-n];
//     
//     fft1(dd, -1);
//     
//     df.resize(2*N);
//     df2.resize(2*N);
//     C one(0, 1);
//     for (int n=0; n<N; n++)
//     {
//       df[n] = -one * dd[n].real() * C(n, 0);
//       df2[n] = df[n] * one * C(n, 0);
//     }
//     df[N] = -one * dd[N].real() * C(0, 0);
//     df2[N] = df[N] * one * C(0, 0); 
//     for (int n=-(N-1); n<=-1; n++)
//     {
//       int i = N+(n+N-1)+1;
//       df[i] = -one * dd[i].real() * C(n, 0);
//       df2[i] = df[i] * one * C(n, 0);
//     }
//     fft1(df, -1);
//     fft1(df2, -1);
//     
//     std::vector<T> result {f[order/2],
//                           -df[order/2].real() / (2*N) * con, 
//                           -df2[order/2].real() / (2*N) * con * con,
//                           0, 0};    
//     
//     mode = central ? 2 : 3;
//     x = central ? sigma : fmean * fmean + sigma;
//     
//     for (size_t n=0; n<order; n++)
//     {
//       if (n != N/2)
//         f[n] = evaluate(z[n]*bma + x);
//     }
// 
//     for (size_t n=0; n<order; n++)
//       dd[n] = f[n];
//         
//     for (size_t n=0; n<order-2; n++)
//       dd[order+n] = f[order-2-n];
//     
//     fft1(dd, -1);
//     
//     for (int n=0; n<N; n++)
//     {
//       df[n] = -one * dd[n].real() * C(n, 0);
//       df2[n] = df[n] * one * C(n, 0);
//     }
//     df[N] = -one * dd[N].real() * C(0, 0);
//     df2[N] = df[N] * one * C(0, 0); 
//     for (int n=-(N-1); n<=-1; n++)
//     {
//       int i = N+(n+N-1)+1;
//       df[i] = -one * dd[i].real() * C(n, 0);
//       df2[i] = df[i] * one * C(n, 0);
//     }
//     fft1(df, -1);
//     fft1(df2, -1);
//     
//     result[3] = -df[order/2].real() / (2*N) * con;
//     result[4] = -df2[order/2].real() / (2*N) * con * con;
//     
//     return result;
//     
    
  }
    
  ComputeTau<T>* ctau=0;
  ComputeTau<float>* ftau=0;
  ComputeTau<double>* dtau=0;
  DCT<T>* dct;
  DST<T>* dst;
  T fmean, sigma;
  float tau0 = -1;
  double tau1;
  bool central = true;
  size_t max_integration_interval;
  size_t min_integration_interval;
  size_t mode;
  size_t max_iters_hp;  
  T x;
  T dx;
  size_t order;
  bool compute_initial_float=true;
};

template <class T>
class DerivativesDP : public DerivativesHP<T>
{
  // Compute a quick initial guess using bracketing float solver
  // Compute refinement using newton with T precision

public:
  
  using DerivativesHP<T>::fmean;
  using DerivativesHP<T>::sigma;
  using DerivativesHP<T>::ftau;
  using DerivativesHP<T>::ctau;
  using DerivativesHP<T>::tau0;
  using DerivativesHP<T>::x;
  using DerivativesHP<T>::mode;
  using DerivativesHP<T>::central;
  using DerivativesHP<T>::order;
  using DerivativesHP<T>::dx;
  using DerivativesHP<T>::dct;
  using DerivativesHP<T>::dst;
  
  DerivativesDP(T fmean, T sigma, T dx, size_t order, bool central)
  : DerivativesHP<T>(fmean, sigma, dx, order, central) 
  {};

  DerivativesDP(T fmean, T sigma, T dx, size_t order, bool central,
                size_t imin, size_t imax)
  : DerivativesHP<T>(fmean, sigma, dx, order, central, imin, imax) 
  {};
  
  void init()
  {
    ctau = new ComputeTau<T>(fmean, sigma, 
                             this->min_integration_interval, 
                             this->max_integration_interval,
                             std::numeric_limits<T>::epsilon(), 
                             2*std::numeric_limits<T>::digits);
    
    ftau = new ComputeTau<float>(float(fmean), float(sigma), 2, 8,
                                 5*std::numeric_limits<float>::epsilon(),
                                 2*std::numeric_limits<float>::digits10);
    
  }
  
  bool set_parameters(T fm, T sm)
  {
    if (!ctau)
      init();

    fmean = fm;
    sigma = sm;
    T seg = sm / fm / (1-fm);
    if (seg > (1-std::numeric_limits<T>::epsilon()) or
        seg < std::numeric_limits<T>::epsilon())
    {
      std::cout << "fmean " << fm << " sigma " << sm << " Seg " << seg << std::endl;
      return false;
    }
    else
    {
      bool ok;
      ok = ctau->set_parameters(fmean, sigma);
      ok = ftau->set_parameters(float(fmean), float(sigma));
      return ok;
    }
  }
  
  T compute_tau(T fmean, T sigma)
  {
    set_parameters(fmean, sigma);
    T result = ctau->compute_newton(T(tau0), 4);
  }
  
  T compute_tau()
  {
    if (!ctau)
      init();

    return ctau->compute_newton(T(tau0), 4);
  }
  
  T evaluate(T const z) 
  {
    T fm0 = fmean;
    T sigma0 = sigma;
    T im0 = fmean * fmean + sigma;
    if (this->mode == 0)
    { // Derivative with respect to fmean (sigma constant)
      fmean = z;
    }
    else if (this->mode == 1)
    { // Derivative with respect to fmean (im constant)
      fmean = z;
      sigma = im0 - fmean * fmean;
    }
    else if(this->mode == 2)
    { // Derivative with respect to sigma (fmean constant)
      sigma = z;
    }
    else if(this->mode == 3)
    {  // Derivative with respect to im (fmean constant)
      sigma = z - fmean * fmean;
    }
    bool ok = set_parameters(fmean, sigma);
    if (not ok)
    {
      return -1;
    }
    else
    {    
      T seg = sigma / fmean / (1-fmean);
      T result;
      if (seg > 0.95)
      {
        result = ctau->compute_bracket(0.0001, 0.1);
      }      
      else
      {
        //if (tau0 < 0)
          
        tau0 = ftau->compute_bracket();    
        
        result = ctau->compute_newton(T(tau0), 4);
      }
      ctau->set_parameters(fm0, sigma0); 
      fmean = fm0;
      sigma = sigma0;
      return result;
    }
  }    
};

#endif
