#include <boost/multiprecision/cpp_dec_float.hpp>
#include "pmf_multiprecision.hpp"

#include <iostream>
#include <iomanip>
using std::setw; 
using std::setprecision; 
using std::scientific;
using std::cout; using std::endl;
using boost::multiprecision::cpp_dec_float_100;  
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float;
using boost::multiprecision::number;

int main()
{  
//   auto_cpu_timer timer;
  typedef number<cpp_dec_float<8> > cpp_dec_float_8;
  typedef number<cpp_dec_float<16> > cpp_dec_float_16;
  typedef number<cpp_dec_float<24> > cpp_dec_float_24;
  typedef number<cpp_dec_float<32> > cpp_dec_float_32;
  typedef number<cpp_dec_float<48> > cpp_dec_float_48;
  typedef number<cpp_dec_float<64> > cpp_dec_float_64;
  
  typedef cpp_dec_float_32 T;
  //typedef double T;
  
  T pi = boost::math::constants::pi<T>();
   T fmean = T(1.0) - T(2) / 100;
//     T fmean = T(2.0) / 100;
  //T fmean = 5.4366172667797141e-02; 
   T sigma = (T(1.0) - T(1) / 1000) * fmean * (T(1.0) - fmean);

   //T sigma = T(2.0) / 100 * fmean * (T(1.0) - fmean);
//   T sigma = 5.4366172667797141e-02* fmean * (T(1.0) - fmean);
  T im = fmean*fmean + sigma;
    
   cout << scientific << setprecision(numeric_limits<T>::digits10);
//    cout << numeric_limits<T>::digits << endl;
//   cout << fmean << " " << sigma << endl;
  T h = T(1.0) / 1000000;  

//   ClenshawCurtisTransformed<T>* CT = new ClenshawCurtisTransformed<T>(std::numeric_limits<T>::epsilon(), fmean, sigma, 2, 12);
//   ClenshawCurtis<T>* CC = new ClenshawCurtis<T>(std::numeric_limits<T>::epsilon(), fmean, sigma, 2, 12);
//   boost::math::tuple<T, T> result0 = CT->integrate_fdf(T(1) / 100);
//   boost::math::tuple<T, T> result1 = CC->integrate_fdf(T(1) / 100);
//   cout << std::get<0>(result0) << " " << std::get<1>(result0) << endl;
//   cout << std::get<0>(result1) << " " << std::get<1>(result1) << endl;
//   cout << CC->integrate(0.25) << endl;
//  cout << result1[0] << " " << result1[1] << " " << result1[2] << endl;
//  T tau = findtau_mp_cc<T>(CT);
//  cout << tau << endl;
//      cout << CT->integrate(tau) << endl;     
     
//     ComputeTau<T>* ctau = new ComputeTau<T>(fmean, sigma, 2, 12);
//     T tau0 = ctau->compute_bracket();
//    T tau1 = ctau->compute_newton(0.99*tau0);
//    cout << tau0 << endl;
//    cout << tau1 << endl;
//    cout << ctau->CC->integrate(tau1) << endl;
//    cout << ctau->CC->integrate(tau0) << endl;
//     boost::math::tuple<T, T> tt = ctau->CC->integrate_fdf(tau1);
//     cout << std::get<0>(tt) << " " << std::get<1>(tt) << endl;


//   DerivativesDP<double>* dder = new DerivativesDP<double>(double(fmean), double(sigma), double(h), 11, true, 2, 8);
//   std::vector<double> result4 = dder->derivatives();
//   cout << result4[0] << " " << result4[3] << " " << result4[4] << endl; 
// 
//   cout << scientific << setprecision(numeric_limits<T>::digits10);
  DerivativesHP<T>* der = new DerivativesHP<T>(fmean, sigma, h, 11, true, 2, 12, 2);
//  der->derivatives();
  der->init();
  std::cout << der->dtau->compute_bracket() << std::endl;
  
  for (size_t i = 5; i<9; i++)
  {
    h = T(1.0) / pow(10, i);  
    der->dx = h;
    std::vector<T> result3 = der->derivatives();
    cout << result3[0] << " " << result3[3] << " " << result3[4] << endl;    
  }
  
  
}
