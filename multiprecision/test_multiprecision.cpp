#include <boost/multiprecision/cpp_dec_float.hpp>
#include "pmf_multiprecision.hpp"

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
  cpu_timer timer2;
  typedef number<cpp_dec_float<8> > cpp_dec_float_8;
  typedef number<cpp_dec_float<16> > cpp_dec_float_16;
  typedef number<cpp_dec_float<24> > cpp_dec_float_24;
  typedef number<cpp_dec_float<32> > cpp_dec_float_32;
  typedef number<cpp_dec_float<48> > cpp_dec_float_48;
  typedef number<cpp_dec_float<64> > cpp_dec_float_64;
  
  typedef cpp_dec_float_32 T;
  //typedef double T;
  
  T pi = boost::math::constants::pi<T>();
   //T fmean = T(1.0) - T(5) / 100;
   T fmean = T(50) / 100;
  //T fmean = 5.4366172667797141e-02; 
   T sigma = (T(1) - T(1) / 10) * fmean * (T(1) - fmean);

   //T sigma = T(2.0) / 100000 * fmean * (T(1.0) - fmean);
//   T sigma = 5.4366172667797141e-02* fmean * (T(1.0) - fmean);
  T im = fmean*fmean + sigma;
   cout << scientific << setprecision(numeric_limits<T>::digits10);
    
//    cout << numeric_limits<T>::digits << endl;
//   cout << fmean << " " << sigma << endl;
  T h = T(1.0) / 1000000; 
  
  ComputeTau<T>* CT = new ComputeTau<T>(fmean, sigma, 2, 10);
  ComputeTau<double>* CD = new ComputeTau<double>(double(fmean), double(sigma), 2, 10);
//   cout << CT->compute_bracket() << endl;
//   cout << CD->compute_bracket() << endl;
//   cout << " " << CT->tau << " " << CT->CC->integrate(CT->tau) << endl;
//   cout << " " << CD->tau << " " << CD->CC->integrate(CD->tau) << endl;
  
  
//    ClenshawCurtisTransformed<T>* CT = new ClenshawCurtisTransformed<T>(std::numeric_limits<T>::epsilon(), fmean, sigma, 2, 12);
// // //    ClenshawCurtis<T>* CC = new ClenshawCurtis<T>(std::numeric_limits<T>::epsilon(), fmean, sigma, 2, 12);    
//    timer2.start();
//    for (size_t n=2; n<10; n++)
//      CT->compute_weights(n);
//    cout << timer2.format() << endl;
//    timer2.start();
//    CT->read_weights();
//    cout << timer2.format() << endl;
//    
   
   
//      std::ifstream myfile32("weights32.dat");
// //     std::ifstream myfile64("weights64.dat");
// //     CT->compute_and_write_weights_to_file(myfile);    
// //    cout << numeric_limits<T>::digits10 << endl;
//      CT->read_weights_from_file(myfile32);
// //     CT->read_weights_from_file(myfile64);
        
//   std::string line;
//   std::ifstream myfile("weights64.dat");
//   double weight;
//   if (myfile.is_open())
//   {
//     while ( std::getline (myfile,line) )
//     {
//       weight = boost::lexical_cast<double>(line);
//       cout << line << '\n';
//       cout << weight << endl;
//     }
//     myfile.close();
//   }
//    boost::math::tuple<T, T> result0;
//    timer2.start();
//    result0 = CT->integrate_fdf(T(1) / 100);
//    cout << timer2.format() << endl;
//    
//    for (size_t i=0; i<10; i++)
//      result0 = CT->integrate_fdf(T(1) / 100);
// //    boost::math::tuple<T, T> result1 = CC->integrate_fdf(T(1) / 100);
//    cout << timer2.format() << endl;
//    cout << std::get<0>(result0) << " " << std::get<1>(result0) << endl;
//    cout << std::get<0>(result1) << " " << std::get<1>(result1) << endl;
//   cout << CC->integrate(0.25) << endl;
//  cout << result1[0] << " " << result1[1] << " " << result1[2] << endl;
//  T tau = findtau_mp_cc<T>(CT);
//  cout << tau << endl;
//      cout << CT->integrate(tau) << endl;
//   cout << std::numeric_limits<float>::epsilon() << endl;
//   cout << std::numeric_limits<double>::epsilon() << endl;
//   cout << std::numeric_limits<T>::epsilon() << endl;   
//     ComputeTau<T>* ctau = new ComputeTau<T>(fmean, sigma, 2, 12);
//     T tau0 = ctau->compute_bracket();
//    T tau1 = ctau->compute_newton(0.99*tau0);
//    cout << tau0 << endl;
//    cout << tau1 << endl;
//    cout << ctau->CC->integrate(tau1) << endl;
//    cout << ctau->CC->integrate(tau0) << endl;
//     boost::math::tuple<T, T> tt = ctau->CC->integrate_fdf(tau1);
//     cout << std::get<0>(tt) << " " << std::get<1>(tt) << endl;
  
  DerivativesDP<double>* dder = new DerivativesDP<double>(double(fmean), double(sigma), double(h), 11, true, 2, 12);
  timer2.start();
  std::vector<double> result4 = dder->derivatives();
  cout << timer2.format() << endl;
  timer2.start();
  result4 = dder->derivatives();
  cout << timer2.format() << endl;
  timer2.start();
  result4 = dder->derivatives();
  cout << timer2.format() << endl;

  timer2.start();
  float t0 = dder->ftau->compute_bracket();
  dder->tau0 = t0;
  double result5 = dder->compute_tau();
  cout << timer2.format() << endl;
  timer2.start();
  result5 = dder->compute_tau();
  cout << timer2.format() << endl;
  timer2.start();
  result5 = dder->compute_tau();
  cout << timer2.format() << endl;

  //   
   cout << result5 << endl;
   cout << result4[0] << " " << result4[1] <<" " << result4[2] << " " << result4[3] << " " << result4[4] << endl;   
// // 
// // //   cout << scientific << setprecision(numeric_limits<T>::digits10);
  DerivativesHP<T>* der = new DerivativesHP<T>(fmean, sigma, h, 11, true, 2, 12, 4);
  timer2.start();
  std::vector<T> result6 = der->derivatives();
  cout << timer2.format() << endl;
  timer2.start();
  result6 = der->derivatives();
  cout << timer2.format() << endl;

  cout << result6[0] << " " << result6[3] << " " << result6[4] << endl; 
// //   der->init();
// //   std::cout << der->ftau->compute_bracket() << std::endl;
// //   std::cout << der->dtau->compute_bracket() << std::endl;
// //   std::cout << der->ctau->compute_bracket() << std::endl;
//   timer2.start();
//   std::vector<T> result3;
//   for (size_t i = 6; i<7; i++)
//   {
//     h = T(1.0) / pow(10, i);  
//     der->dx = h;
//     result3 = der->derivatives();
//     cout << result3[0] << " " << result3[1] << " " << result3[2] << " " << result3[3] << " " << result3[4] << endl;    
//   }
//   cout << timer2.format() << endl;
//   
// //   cout << policies::get_max_root_iterations<Policy>();
//   cout << std::numeric_limits<float>::digits << endl;
//   cout << std::numeric_limits<float>::digits10 << endl;
//   cout << std::numeric_limits<double>::digits << endl;
//   cout << std::numeric_limits<double>::digits10 << endl;
//   cout << std::numeric_limits<T>::digits << endl;
//   cout << std::numeric_limits<T>::digits10 << endl;
// 
  
}
