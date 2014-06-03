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
// First added:  2012-05-14
// Last changed: 2012-05-14

#include "PMF.h"
#include <time.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "pmf_multiprecision.hpp"
#include <boost/fusion/include/std_pair.hpp>

using namespace pmfpack;
using namespace std;
using boost::multiprecision::number;
using boost::multiprecision::cpp_dec_float;
using std::setprecision; 
using std::scientific;
using std::cout; using std::endl;

void printvec(double* x, int N, string s)
{
  std::cout << s << std::endl;
  for (int i=0; i<N; i++)
  {
    std::cout << " " << i << " " << x[i] << std::endl; 
  }
}

int main()
{
    typedef number<cpp_dec_float<32> > cpp_dec_float_32;
    typedef number<cpp_dec_float<48> > cpp_dec_float_48;
    
   typedef cpp_dec_float_32 T;
//    typedef double T;
    
//     double fmean, sigma, intensity_of_segregation;
    
    // Use some numbers from a 1D simulation of a mixing layer
//     fmean = 0.02;
    T fmean = T(2.0) / 1000;
//     T sigma = (T(1.0) - T(6.0) / 10) * fmean * (T(1.0) - fmean);
    T sigma = T(1.0) / 100 * fmean * (T(1.0) - fmean);
//     intensity_of_segregation = 0.4;
//     sigma = intensity_of_segregation * fmean * (1 - fmean);
//     sigma = 0.2559307963e-1; 
    
    T h = T(1.0) / 100000;
//    DerivativesHP<T>* der = new DerivativesHP<T>(T(fmean), T(sigma), h, 11, true, 2, 9);
    DerivativesDP<T>* der = new DerivativesDP<T>(T(fmean), T(sigma), h, 11, true, 2, 9);
    std::vector<T> result3 = der->derivatives();
    cout << scientific << setprecision(numeric_limits<T>::digits10);
    cout << result3[0] << " " << result3[1] << " " << result3[3] << " " << result3[2] << " " << result3[4] << endl;
    
    // Create a new PMF object 
    PMF *pmf = new PMF(double(fmean), double(sigma), 1, 1, 0, 1);
    
    // Specify mean flow parameters from the simulation
    pmf->set_fmean_gradient(0.2829907608e+2, -0.5249074097e+3, 0);
    pmf->set_sigma_gradient(0.2471093655e+1, -0.5286500931e+2, 0);
    pmf->chi = 0.2305485382e+3;
    pmf->DT = 0.1144416248e-3;
    pmf->central = true;
    // Compute tau and derivatives of tau
//      pmf->roots->froot->realloc(2);
    pmf->roots->froot->error_message_off();
    pmf->verbose = 0;
    clock_t start = clock();
    for (int i=0; i<100; i++){
      pmf->tau = 0.1;
      pmf->compute(0, true);
    }
    double duration = (double)(clock() - start) / CLOCKS_PER_SEC;
    cout << setprecision(numeric_limits<double>::digits10);
    std::cout << "Computed tau  " << pmf->tau << " time " << duration/100 << std::endl;
    std::cout << pmf->tau << " " << pmf->dtaudf << " " << pmf->dtauds <<  " " << pmf->d2taudfdf << " " << pmf->d2taudsds << std::endl;

//     start = clock();
//     for (int i=0; i<100; i++){
//       der->compute_tau(fmean, sigma);
//     }
//     duration = (double)(clock() - start) / CLOCKS_PER_SEC;
//     std::cout << "time " << duration/100 << std::endl;
    // Test lookuptable
     HDF5Lookup *lookuptable = new HDF5Lookup(pmf->derivator);
     lookuptable->generate_table(50, 50, "HDF5_50_24_table.h5");
     //lookuptable->read_table("HDF5_10_table.h5");
     //lookuptable->set_order(6);
//    GSLLookup *lookuptable = new GSLLookup(pmf->derivator);
//    lookuptable->generate_table(10, 10, "GSL_10_table.dat");
//     
//    lookuptable->read_table("GSL_test_table.dat");

//     std::cout << "Look up tau etc" << std::endl;
//     start = clock();
//     for (int i=0; i<1; i++){
//       lookuptable->compute(0, true);  // Look up tau etc
// //       lookuptable->compute(0, true, true);  // With root-polishing
//     }
//     duration = (double)(clock() - start) / CLOCKS_PER_SEC;
//     std::cout << "Looked up tau  " << pmf->tau << " time " << duration/100000. << std::endl;
//         
//     T h = T(1.0) / pow(10, 5);
//     std::vector<T> ctau1 = chebfft(T(sigma), h, 11, Boost_dx_Function<T>(CC, 2, findtau_mp_cc<T>));
//     std::vector<T> ctau2 = chebfft(T(fmean), h, 11, Boost_dx_Function<T>(CC, 0, findtau_mp_cc<T>));
//     
//     std::cout << ctau1[0] << " " << ctau2[1] << " " << ctau1[1] << " " << ctau2[2] << " " << ctau1[2] << std::endl;
//     
//     std::cout << pmf->tau << " " << pmf->dtaudf << " " << pmf->dtauds <<  " " << pmf->d2taudfdf << " " << pmf->d2taudsds << std::endl;

//     // Compute all possible conditional models
//     int N = 86;
//     double *eta = (double *) malloc(N * sizeof(double));
//     double *p = (double *) malloc(N * sizeof(double));
//     double *cv = (double *) malloc(3 * N * sizeof(double));
//     for (int i=1; i<N; i++)
//       eta[i] = (i + 1) * 1. / (N+1);    
// 
//     printvec(eta, N, "eta");
//     pmf->counterflow(eta, N, p, N);
//     printvec(p, N, " counterflow");
//     pmf->CSD(eta, N, p, N);
//     printvec(p, N, " CSD");
//     pmf->PDF(eta, N, p, N);
//     printvec(p, N, " PDF");
//     pmf->Laplace(eta, N, p, N);
//     printvec(p, N, " Laplace");
//     pmf->CV(eta, N, cv, 3, N);
//     printvec(cv, N, " CV[0]");
//     std::cout << " counterflow " << pmf->counterflow(0.2) << std::endl; 
//     std::cout << " CSD " << pmf->CSD(0.2) << std::endl;
//     std::cout << " PDF " << pmf->PDF(0.2) << std::endl;
//     std::cout << " Laplace " << pmf->Laplace(0.2) << std::endl;
//     pmf->CV(0.2, cv, 3);
//     std::cout << " CV[0] " << cv[0] << std::endl;
//     start = clock();
//     for (int i=0; i<100000; i++) pmf->CSD_verify(eta, N, p, N);
//     duration = (double)(clock() - start) / CLOCKS_PER_SEC;
//     std::cout << "Time CSD2 " << duration / 100000. << " value " << p[4] << endl;
//     start = clock();
//     for (int i=0; i<100000; i++) pmf->CSD(eta, N, p, N);
//     duration = (double)(clock() - start) / CLOCKS_PER_SEC;
//     std::cout << "Time CSD  " << duration / 100000. << " value " << p[4] << endl;
}
