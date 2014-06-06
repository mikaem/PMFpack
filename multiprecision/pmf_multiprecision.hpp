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

using std::numeric_limits; 
using boost::timer::auto_cpu_timer;

inline bool is_pow_2(size_t n) {
  return n & (n - 1) == 0;
}

inline size_t next_pow_2(size_t x) {
  size_t h = 1;
  while(h < x) h *= 2;
  return h;
}

// exp(i 2 pi m / n)
template<class N>
static std::complex<N> exp(int m, int n) {
  static const N pi = boost::math::constants::pi<N>();
  return std::polar<N>(1, pi * 2 * m / n);
}

// bit-reverse transpose the array
template<class T>
static void bitrev(std::vector<T> &a) {
  const size_t n = a.size();
  for(size_t i = 0, j = 0 ; i < n - 1 ; i++) {
    if(i < j) std::swap(a[i], a[j]);
    // bit reversed counter
    size_t m = n;
    do {
        m >>= 1;
        j ^= m;
    } while(!(j & m));
  }
}

/// compute a power-of-2 FFT in-place.
template<class N>
static void fft0(std::vector<std::complex<N>> &a, int sign) {
  const size_t n = a.size();
  bitrev(a);
  for(size_t i = 1 ; i < n ; i *= 2)
  for(size_t j = 0 ; j < i ; j++) {
    const auto w = exp<N>(sign * j, 2 * i);
      for(size_t k = j ; k < n ; k += 2 * i) {
        const auto v = a[k];
        const auto x = a[k + i] * w;
        a[k] = v + x;
        a[k + i] = v - x;
    }
  }
}

/// compute a Bluestein Chirp-Z transform in-place.
template<class N>
void bluestein(std::vector<std::complex<N>> &a) {
  const size_t n = a.size();
  const size_t nb = next_pow_2(2 * n);
  const std::complex<N> zero(0, 0);
  std::vector<std::complex<N>> w(n), y(nb, zero), b(nb, zero);
  for(size_t k = 0 ; k < n ; k++)
    w[k] = exp<N>(((unsigned long long)k * k) % (2 * n), 2 * n);
  for(size_t i = 0 ; i < n ; i++) y[i] = w[i];
  for(size_t i = 1 ; i < n ; i++) y[nb - i] = w[i];
  fft0(y, -1);
  for(size_t i = 0 ; i < n ; i++) b[i] = conj(w[i]) * a[i];
  // scaled convolution b * y
  fft0(b, -1);
  for(size_t i = 0 ; i < nb ; i++) b[i] *= y[i];
  fft0(b, 1);
  for(size_t i = 0 ; i < n ; i++) a[i] = conj(w[i]) * b[i] / N(nb);
}

/// Swap real and imaginary parts of all elements.
template<class N>
static void swapri(std::vector<std::complex<N>> &a) {
  for(auto &x : a) x = std::complex<N>(std::imag(x), std::real(x));
}

/// Compute an arbitrary-sized FFT in-place.
template<class N>
static void fft1(std::vector<std::complex<N>> &a, int sign) {
  if(is_pow_2(a.size()))
    fft0(a, sign);
  else {
    if(sign == 1) swapri(a);
    bluestein(a);
    if(sign == 1) swapri(a);
  }
}

/// Compute an arbitrary-sized FFT out-of-place.
template<class N, class T>
static std::vector<std::complex<N>> simple_fft(const std::vector<T> &a, int sign = -1) {
    const int n = a.size();
    std::vector<std::complex<N>> b(a.begin(), a.end());
    fft1(b, sign);
    if(sign == 1) for(auto &x : b) x /= n;
    return b;
}

template <class T>
class DCT
{ // Discrete cosine transform type I
public:
  DCT(size_t N)
  : N(N)
  {
    init();
  }
  
  void init()
  {
    static const T pi = boost::math::constants::pi<T>();
    for (size_t k=0; k<N; k++){
      for (size_t n=1; n<N-1; n++){
        size_t m = (n*k) % ((N-1)/2);
        typename std::map<size_t, T>::iterator it = cosmap.find(m);
        if (it == cosmap.end())
          cosmap[m] = cos(pi*(T(m) / (N-1)));
      }
    }
  }
  
  std::vector<T> dct_odd(const std::vector<T> &x)
  {
    static const T half = boost::math::constants::half<T>();
    assert(N % 2 == 1);
    assert(N == x.size());
    
    std::vector<T> xk(N);  
    int sig = -1;
    for (size_t k=0; k<N; k++){
      sig *= -1;
      xk[k] = half*(x[0]+sig*x[N-1]);
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
    static const T half = boost::math::constants::half<T>();
    assert(N % 2 == 1);
    assert(N == x.size());
    int sig = pow(-1, k);
    T xk = half*(x[0]+sig*x[N-1]);
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
        typename std::map<size_t, T>::iterator it = sinmap.find(m);
        if (it == sinmap.end())
          sinmap[m] = sin(pi*(T(m) / (N+1)));
      }
    }
  }  
  
  std::vector<T> dst_odd(const std::vector<T> &x)
  {
  //   auto_cpu_timer timer;
    static const T half = boost::math::constants::half<T>();
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
    static const T half = boost::math::constants::half<T>();
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
    for (size_t i = interval_min; i < interval_max; i++)
      compute_weights(i);    
  }

  ClenshawCurtisBase(T epsilon, T fmean, T sigma, size_t imin, size_t imax) 
  : epsilon(epsilon), fmean(fmean), sigma(sigma), interval_min(imin), interval_max(imax)
  {
    im = fmean * fmean + sigma;
    for (size_t i = interval_min; i < interval_max; i++)
      compute_weights(i);
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
        
        std::vector<T> xx = this->x[i];
        for (size_t n=0; n<Z; n++)
          y[n] = xx[n]*xmhalf + xmsum;
        
        if (i == this->interval_min)
        {
          for (size_t n=0; n<Z; n++)
            fx[n] = function(this->y[n], tau, alfa);
        }
        else
        { // Reuse function values from coarser grid
          for (size_t n=1; n<Z; n=n+2)
            fx[n] = function(y[n], tau, alfa);
          for (size_t n=0; n<Z; n=n+2)
            fx[n] = previous_fvals[n/2];
        }
        // Store function evaluations for even finer grid
        previous_fvals.resize(Z);
        for (size_t n=0; n<Z; n++)
          previous_fvals[n] = fx[n];
        
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << this->fx[n] << std::endl;
        
        new_integral = 0;
        std::vector<T> ww = weights[i];
        for (size_t n=0; n<Z; n++)
          new_integral = new_integral + ww[n]*fx[n];   
        new_integral = new_integral*xmhalf;
        
//           std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
//          std::cout << i << " " << fabs(new_integral - prev_integral) << std::endl;
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
    size_t Z = pow(2, i) + 1; 
    static const T pi = boost::math::constants::pi<T>();
    size_t M = Z-1;
    std::vector<T> z(Z);
    std::vector<T> w(Z);
    
    for (size_t n=0; n<Z; n++)
    {   
      int i = n*2-M;
      z[n] = sin(i*pi/(2*M));
    }
    x[i] = z;
    
    c.resize(Z);
    for (size_t n=0; n<Z; n=n+2)
      c[n] = T(2) / (T(1) - T(n*n));
    
    DCT<T>* dct = new DCT<T>(Z);
    std::vector<T> ff = dct->dct_odd(c);
    std::vector<T> ww(Z);
    w[0] = ff[0] / M;
    w[Z-1] = ff[Z-1] / M;
    for (size_t n=1; n<Z-1; n++)
      w[n] = 2*ff[n] / M;
        
    weights[i] = w;
    
//     std::cout << i << std::endl;
//     for (size_t n=0; n<w.size(); n++)
//     {
//       std::cout << " w "<< n << " " << w[n] << std::endl;
//     }
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
  std::map<size_t, std::vector<T> > weights; 
  std::vector<T> fx, y, c;
};

template <class T> 
class ClenshawCurtisTransformed: public ClenshawCurtisBase<T>
{
public:
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
    T alfa = sqrt(T(2.0)) * boost::math::erfc_inv(2*this->fmean);
    // Compute constant divisor in integral
    T div = 4*sqrt(pi*(2-4*tau));
    T xmin;
    T xmax;
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
//         std::cout << xmin << " " << xmax <<  " " << a << " " << b << std::endl;      
      func = 0;
      integral = this->integrate_domain(xmin, xmax, tau, alfa); 
    }
    else
    {    
      xmin = T(1.0) / (1+2*sqrt(tau)*small);
      func = 2;
      integral = this->integrate_domain(xmin, T(1.0), tau, alfa); 
      
      xmax = sqrt(2*(1-2*tau)*small2);
      xmax = (xmax-alfa) / (1+xmax-alfa);      
//         std::cout << xmin << " " << xmax << std::endl;
      func = 1;
      integral  += this->integrate_domain(T(0.0), xmax, tau, alfa); 
    }
    
//       std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
//       std::cout <<  tau << " " << T(1.0) - integral / this->im / div << std::endl;

    return T(1.0) - integral / this->im / div;    
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
      T df = X * (alfa - phi / (T(1)-T(2)*tau)) * exp(-tmp*tmp-phi*phi/(T(2)-T(4)*tau));
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
    T alfa = sqrt(T(2.0)) * boost::math::erfc_inv(2*this->fmean);
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
      xmin = T(1.0) / (1+2*sqrt(tau)*small);
      func = 2;
      integrals = this->integrate_domain_fdf(xmin, T(1.0), tau, alfa); 
            
      xmax = sqrt(2*(1-2*tau)*small2);
      xmax = (xmax-alfa) / (1+xmax-alfa);      
      
      func = 1;
      std::vector<T> fdintegrals  = this->integrate_domain_fdf(T(0.0), xmax, tau, alfa); 
      integrals[0] += fdintegrals[0];
      integrals[1] += fdintegrals[1];
    }
    
    // Scale the integrals appropriately
    integrals[1] = integrals[1] / tau / sqrt(tau) / root_pi;
    integrals[1] = -integrals[1] / this->im / div;
    
    integrals[0] = T(1.0) - integrals[0] / this->im / div;
    
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
    
    for (size_t i = this->interval_min; i < this->interval_max; i++)
    {
      prev_integral = new_integral;

      size_t Z = pow(2, i) + 1;
      
      fdf.resize(Z);
      this->fx.resize(Z);
      this->y.resize(Z);
      
      std::vector<T> xx = this->x[i];
      for (size_t n=0; n<Z; n++)
        this->y[n] = xx[n]*xmhalf + xmsum;
      
      if (i == this->interval_min)
      {
        for (size_t n=0; n<Z; n++)
        {
          std::vector<T> vals = function_fdf(this->y[n], tau, alfa);
          this->fx[n] = vals[0];
          fdf[n] = vals[1];
        }
      }
      else
      { // Reuse function values from coarser grid
        for (size_t n=1; n<Z; n=n+2)
        {
          std::vector<T> vals = function_fdf(this->y[n], tau, alfa);
          this->fx[n] = vals[0];
          fdf[n] = vals[1];
        }
        for (size_t n=0; n<Z; n=n+2)
        {
          this->fx[n] = previous_fvals[n/2];
          fdf[n] = previous_fdfvals[n/2];
        }
      }
      // Store function evaluations for even finer grid
      previous_fvals.resize(Z);
      previous_fdfvals.resize(Z);
      for (size_t n=0; n<Z; n++)
      {
        previous_fvals[n] = this->fx[n];
        previous_fdfvals[n] = fdf[n];
      }
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << this->fx[n] << " " << fdf[n] << std::endl;
      
      new_integral = T(0);
      std::vector<T> ww = this->weights[i];
      for (size_t n=0; n<Z; n++)
        new_integral = new_integral + ww[n]*this->fx[n];   
      new_integral = new_integral*xmhalf;

      fdf_integral = T(0);
      for (size_t n=0; n<Z; n++)
        fdf_integral = fdf_integral + ww[n]*fdf[n];   
      fdf_integral = fdf_integral*xmhalf;
      
//           std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
//          std::cout << i << " " << fabs(new_integral - prev_integral) << std::endl;
      if (fabs(new_integral - prev_integral) < this->epsilon)
        break;        
    }      
    std::vector<T> result = {new_integral, fdf_integral};
    return result;
  }
  
  std::vector<T> fdf;
  T integral;
  T small = boost::math::erfc_inv(this->epsilon);
//   T small2 = -log(this->epsilon/10000000);
  T small2 = 6*std::numeric_limits<T>::digits;
  size_t func;
};

template <class T> 
class ClenshawCurtis: public ClenshawCurtisBase<T>
{
public:
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
  
  T integrate(T tau)
  {
    static const T pi = boost::math::constants::pi<T>();
    T alfa = sqrt(T(2.0)) * boost::math::erfc_inv(2*this->fmean);
    T integral = 0.0;
    T integral_dom = 1.0;
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
    T xmax = 0.0;
    size_t i = 1;
    
    T dphi = (phi_X_99-phi_X_01) / 2;
    dphi = dphi < T(2.0) / 10 ? dphi : T(2.0) / 10;
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
    integral_dom = 1.0;
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
    
    return T(1.0) - integral / this->im;    
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
        
        this->fx.resize(Z);
        this->y.resize(Z);
        
        std::vector<T> xx = this->x[i];
        for (size_t n=0; n<Z; n++)
          this->y[n] = xx[n]*xmhalf + xmsum;
        
        if (i == this->interval_min)
        {
          if (right and f_upper > this->epsilon)
          {
            this->fx[0] = f_upper;
            for (size_t n=1; n<Z; n++)
              this->fx[n] = function(this->y[n], tau, alfa);
          }
          else if (not right and f_lower > this->epsilon)
          {
            this->fx[Z-1] = f_lower;
            for (size_t n=0; n<Z-1; n++)
              this->fx[n] = function(this->y[n], tau, alfa);
          }
          else
          {           
            for (size_t n=0; n<Z; n++)
              this->fx[n] = function(this->y[n], tau, alfa);
          }
        }
        else
        { // Reuse function values from coarser grid
          for (size_t n=1; n<Z; n=n+2)
            this->fx[n] = function(this->y[n], tau, alfa);
          for (size_t n=0; n<Z; n=n+2)
            this->fx[n] = previous_fvals[n/2];
        }
        // Store function evaluations for even finer grid
        previous_fvals.resize(Z);
        for (size_t n=0; n<Z; n++)
          previous_fvals[n] = this->fx[n];
        
//         for (size_t n=0; n<Z; n++)
//           std::cout << " i " << n << " " << fx[n] << std::endl;
        
        // Store edge values for next interval
        f_upper = this->fx[Z-1];
        f_lower = this->fx[0];
         
        new_integral = 0;
        std::vector<T> ww = this->weights[i];
        for (size_t n=0; n<Z; n++)
          new_integral = new_integral + ww[n]*this->fx[n];   
        new_integral = new_integral*xmhalf;
        
//          std::cout << i << " " << xmin << " " << xmax << " " << new_integral <<  " "<< prev_integral << " " << fabs(new_integral - prev_integral) << std::endl;
        if (fabs(new_integral - prev_integral) < this->epsilon)
          break;        
    }        
    return new_integral;
  }
  
    
  boost::math::normal_distribution<T> s;
  T inv_eps = -quantile(s, numeric_limits<T>::epsilon()/1000);
  T inv_99  = -quantile(s, T(1.0) / 100);
  T inv_01  = -quantile(complement(s, T(1.0) / 100));
  bool right;
  T f_upper, f_lower;
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
//        std::cout << std::scientific << std::setprecision(std::numeric_limits<T>::digits10);
//        std::cout << " " << std::get<0>(tmp) << " " << std::get<1>(tmp) << std::endl;      
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
    digits = std::numeric_limits<T>::digits;
    tol = new tolerance<T>(epsilon);
    init();
  }

  ComputeTau(T fm, T sm, size_t imin, size_t imax)
  :  fmean(fm), sigma(sm), min_integration_interval(imin), 
     max_integration_interval(imax)
  {
    epsilon = std::numeric_limits<T>::epsilon();
    digits = std::numeric_limits<T>::digits;
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
  
  void set_parameters(T fm, T sm)
  {
    fmean = fm;
    sigma = sm;
    CC->set_parameters(fm, sm);
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
  
  T compute_bracket()
  {  
    it = 100;
    T lower = T(1.0) / 100;
    T upper = T(5.0) / 10 - lower;    

    T fa = f->operator()(lower);
    T fb = f->operator()(upper);
    while (fa * fb > 0)
    {
      lower = lower / 10;
      upper = T(5.0) / 10 - lower; 
      fa = f->operator()(lower);
      fb = f->operator()(upper);
    }
    std::pair<T, T> found = boost::math::tools::toms748_solve(*f, lower, upper, fa, fb, *tol, it);
    if (it == 100) std::cout << " Root not found Boost toms748 Fsolver" << std::endl;
    T tau = (found.first + found.second) / 2;
    return tau;
  }
    
  ClenshawCurtisTransformed<T>* CC;
  BoostFunction<T>* f;
  Boost_df_Function<T>* fdf;
  boost::uintmax_t it;  
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
    ctau->set_parameters(fm, sm);
    ftau->set_parameters(M(fm), M(sm));
    M tau0 = ftau->compute_bracket();      
    return ctau->compute_newton(T(tau0));
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
                              2*std::numeric_limits<T>::digits);
    
    dtau = new ComputeTau<double>(double(fmean), double(sigma), 2, 8);
    ftau = new ComputeTau<float>(float(fmean), float(sigma), 2, 7);
    tau0 = ftau->compute_bracket();
  }
  
  virtual void set_parameters(T fm, T sm)
  {
    if (tau0 < 0)
      init();

    fmean = fm;
    sigma = sm;
    ctau->set_parameters(fmean, sigma);
    ftau->set_parameters(float(fmean), float(sigma));
    dtau->set_parameters(double(fmean), double(sigma));
    tau0 = ftau->compute_bracket(); // Initial guess for double    
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
    
    dtau->set_parameters(double(fmean), double(sigma)); 
    double tau1 = dtau->compute_newton(double(tau0), 3); 
    
    ctau->set_parameters(fmean, sigma);
    T result = ctau->compute_newton(T(tau1), max_iters_hp);
    ctau->set_parameters(fm0, sigma0);
    fmean = fm0;
    sigma = sigma0;
    
    return result;
  };
  
  virtual T compute_tau(T fmean, T sigma)
  {
    set_parameters(fmean, sigma);
    double tau1 = dtau->compute_newton(double(tau0), 3);     
    return ctau->compute_newton(T(tau1), max_iters_hp);
  }
  
  virtual T compute_tau()
  {
    double tau1 = dtau->compute_newton(double(tau0), 3);     
    return ctau->compute_newton(T(tau1), max_iters_hp);
  }
    
  virtual std::vector<T> derivatives()
  {
//     auto_cpu_timer timer;  
    if (tau0 < 0)
      init();

    const T half = boost::math::constants::half<T>();
    const T pi = boost::math::constants::pi<T>();
    T bma = dx * half;
    T con = 1 / bma;
    
    x = fmean;
    mode = central ? 0 : 1;

    int N = order - 1;
    std::vector<T> z(order);
    std::vector<T> f(order);
    for (size_t n=0; n<order; n++)
    {
      z[n] = cos(n*pi/N);
      f[n] = evaluate(z[n]*bma + x);
    }

    std::vector<T> dcf = dct->dct_odd(f);
    for (size_t i=0; i<dcf.size() ; i++)
      dcf[i] = dcf[i] * 2;
    
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
  }
    
  ComputeTau<T>* ctau;
  ComputeTau<float>* ftau;
  ComputeTau<double>* dtau;
  DCT<T>* dct;
  DST<T>* dst;
  T fmean, sigma;
  float tau0 = -1;
  bool central = true;
  bool full_precision = true;
  size_t max_integration_interval;
  size_t min_integration_interval;
  size_t mode;
  size_t max_iters_hp;  
  T x;
  T dx;
  size_t order;
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
    
    ftau = new ComputeTau<float>(float(fmean), float(sigma), 2, 7);
    tau0 = ftau->compute_bracket();
  }
  
  void set_parameters(T fm, T sm)
  {
    fmean = fm;
    sigma = sm;
    ctau->set_parameters(fmean, sigma);
    ftau->set_parameters(float(fmean), float(sigma));
    tau0 = ftau->compute_bracket(); // Initial guess for double
  }
  
  T compute_tau(T fmean, T sigma)
  {
    set_parameters(fmean, sigma);
    T result = ctau->compute_newton(T(tau0), 4);
  }
  
  T compute_tau()
  {
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
    
    ctau->set_parameters(fmean, sigma); 
    T result = ctau->compute_newton(T(tau0), 4);
    ctau->set_parameters(fm0, sigma0); 
    fmean = fm0;
    sigma = sigma0;
    return result;
  };
    
};

#endif
