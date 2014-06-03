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
      
//   T fmean = T(1.0) - T(2.0) / 100;
//     T fmean = T(2.0) / 100;
    T fmean = T(1.0) - 0.024891817906210312;
  //T fmean = 5.4366172667797141e-02; 
//    T sigma = (T(1.0) - T(99) / 100) * fmean * (T(1.0) - fmean);
   T sigma = 0.024891817906210312* fmean * (T(1.0) - fmean);

   //T sigma = T(2.0) / 100 * fmean * (T(1.0) - fmean);
//   T sigma = 5.4366172667797141e-02* fmean * (T(1.0) - fmean);
  T im = fmean*fmean + sigma;
    
   cout << scientific << setprecision(numeric_limits<T>::digits10);
//    cout << numeric_limits<T>::digits << endl;
//   cout << fmean << " " << sigma << endl;
  T h = T(1.0) / 100000;  

//  ClenshawCurtisTransformed<T>* CT = new ClenshawCurtisTransformed<T>(std::numeric_limits<T>::epsilon(), fmean, sigma, 2, 10);
//  std::vector<T> result1 = chebfft(sigma, h, 11, Boost_dx_Function<T>(CT, 2, findtau_mp_cc<T>));
//  cout << result1[0] << " " << result1[1] << " " << result1[2] << endl;
//  T tau = findtau_mp_cc<T>(CT);
//  cout << tau << endl;
//      cout << CT->integrate(tau) << endl;     
     
//    ComputeTau<T>* ctau = new ComputeTau<T>(fmean, sigma, 2, 10);
//    T tau0 = ctau->compute_bracket();
//    T tau1 = ctau->compute_newton(0.99*tau0);
//    cout << tau0 << endl;
//    cout << tau1 << endl;
//    cout << ctau->CC->integrate(tau1) << endl;
//    cout << ctau->CC->integrate(tau0) << endl;
//     boost::math::tuple<T, T> tt = ctau->CC->integrate_fdf(tau1);
//     cout << std::get<0>(tt) << " " << std::get<1>(tt) << endl;

  DerivativesHP<T>* der = new DerivativesHP<T>(T(0.5), T(0.1), h, 11, true, 2, 9);
  der->set_parameters(fmean, sigma);
  std::vector<T> result3 = der->derivatives();
  cout << result3[0] << " " << result3[3] << " " << result3[4] << endl;

  cout << scientific << setprecision(numeric_limits<double>::digits10);
  DerivativesDP<double>* dder = new DerivativesDP<double>(double(fmean), double(sigma), double(h), 11, true, 2, 8);
  std::vector<double> result4 = dder->derivatives();
  cout << result4[0] << " " << result4[3] << " " << result4[4] << endl; 
  
/*
  T var = sigma;
  const T dx11 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 4)*var;
  const T dx21 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 3)*var;
  const T dx31 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 2)*var;  
  const T dx1 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 8)*var;
  const T dx2 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 6)*var;
  const T dx3 = pow(std::numeric_limits<T>::epsilon(), T(1.0) / 4)*var;*/
  
//   for (int i=2; i<=12; i=i+1)
//   {
//     T h = T(1.0) / pow(10, i) * im;
//     std::pair<T, T> result = derivative(im, h, Boost_dx_Function<T>(CC, 3, findtau_mp_cc<T>));
//     cout << h << " " <<  result.first << " " << result.second << endl;
//   }
//   cout << endl;
//   for (int i=3; i<=12; i=i+1)
//   {
//     T h = T(1.0) / pow(10, i);
//     std::pair<T, T> result = chebfft(sigma, h, 11, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//     //std::pair<T, T> result = chebfft(fmean, h, 11, Boost_dx_Function<T>(CC, 0, findtau_mp_cc<T>));
//     cout << h << " " <<  result.first << " " << result.second << endl;
//   }
//   cout << CC->tau << endl;
  
//   T h = T(1.0) / pow(10, 8) * fmean;
//   std::pair<T, T> result = chebfft(fmean, h, 11, Boost_dx_Function<T>(CC, 0, findtau_mp_cc<T>));
//   cout << h << " " <<  result.first << " " << result.second << endl;

//   cout << dx11 << endl;
//   cout << dx21 << endl;
//   cout << dx31 << endl;
//   cout << "Cheb " << endl;
//   cout << dx1 << endl;
//   cout << dx2 << endl;
//   cout << dx3 << endl;
// 
//   std::pair<T, T> dtauds1 = derivative(sigma, dx11, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> dtauds2 = derivative(sigma, dx21, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> dtauds3 = derivative(sigma, dx31, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   cout << " Cheb 1 " << endl;
   
//   std::pair<T, T> cdtauds21 = chebfft(sigma, dx1, 5, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds31 = chebfft(sigma, dx1, 7, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds41 = chebfft(sigma, dx1, 9, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   cout << " Cheb 2 " << endl;
//   std::pair<T, T> cdtauds12 = chebfft(sigma, dx2, 3, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds22 = chebfft(sigma, dx2, 5, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds32 = chebfft(sigma, dx2, 7, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds42 = chebfft(sigma, dx2, 9, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   cout << " Cheb 3 " << endl;
//   std::pair<T, T> cdtauds13 = chebfft(sigma, dx3, 3, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds23 = chebfft(sigma, dx3, 5, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds33 = chebfft(sigma, dx3, 7, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//   std::pair<T, T> cdtauds43 = chebfft(sigma, dx3, 9, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));  
// 
//   cout << "First derivative " << endl;
//     cout << dtauds1 << endl;
//     cout << dtauds2 << endl;
//     cout << dtauds3 << endl;
// 
//   cout << "Cheb 1 first derivative " << endl;
//     cout << cdtauds11 << endl;
//     cout << cdtauds21 << endl;
//     cout << cdtauds31 << endl;
//     cout << cdtauds41 << endl;
//   cout << "Cheb 2 first derivative " << endl;
//     cout << cdtauds12 << endl;
//     cout << cdtauds22 << endl;
//     cout << cdtauds32 << endl;
//     cout << cdtauds42 << endl;
//   cout << "Cheb 3 first derivative " << endl;
//     cout << cdtauds13 << endl;
//     cout << cdtauds23 << endl;
//     cout << cdtauds33 << endl;
//     cout << cdtauds43 << endl;
   
// const T dxx = pow(std::numeric_limits<T>::epsilon(), T(1.0)/3)*boost::math::constants::pi<T>() / 3;
// cout << dxx << endl;
// 
// const T d_mp = derivative( 
//       T(boost::math::constants::pi<T>() / 3),
//       dxx,
//       [](const T& x) -> T
//       {
//          return sin(x);
//       }
//       );
// // const T d_mp2 = chebfft(T(boost::math::constants::pi<T>() / 3),
// //       T(1.0E-5), 3,
// //       [](const T& x) -> T
// //       {return sin(x);});
// // 
// // const T d_mp3 = chebfft(T(boost::math::constants::pi<T>() / 3),
// //       T(1.0E-5), 5,
// //       [](const T& x) -> T
// //       {return sin(x);});
// // const T d_mp4 = chebfft(T(boost::math::constants::pi<T>() / 3),
// //       T(1.0E-5), 7,
// //       [](const T& x) -> T
// //       {return sin(x);});
// // const T d_mp5 = chebfft(T(boost::math::constants::pi<T>() / 3),
// //       T(1.0E-5), 9,
// //       [](const T& x) -> T
// //       {return sin(x);});
// // 
// 
//    cout << -sin(boost::math::constants::pi<T>() / 3) << endl;
//    cout << d_mp << endl;
// //   cout << d_mp2 << endl;
// //   cout << d_mp3 << endl;
// //   cout << d_mp4 << endl;
// //   cout << d_mp5 << endl;

}
