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

using namespace pmfpack;
using namespace std;

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
    double fmean, sigma, intensity_of_segregation;
    
    // Use some numbers from a 1D simulation of a mixing layer
    fmean = 0.70837758;
    intensity_of_segregation = 0.194360748; 
    sigma = intensity_of_segregation * fmean * (1 - fmean);
    
    // Create a new PMF object 
    PMF *pmf = new PMF(fmean, sigma, 0, 0, 0, 2);
    
    // Specify mean flow parameters from the simulation
    pmf->set_fmean_gradient(-0.3496497, 0, 0);
    pmf->set_sigma_gradient(-0.4294126, 0, 0);
    pmf->chi = 1.003770;
    pmf->DT = 1.;

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
//    std::cout.precision(8);
    std::cout << "Computed tau  " << pmf->tau << " time " << duration/100 << std::endl;
    std::cout << pmf->tau << " " << pmf->dtaudf << " " << pmf->dtauds << std::endl;
    
    // Test lookuptable
    GSLLookup *lookuptable = new GSLLookup(pmf->derivator);
    lookuptable->generate_table(10, 10, "GSL_test_table.dat");
    //lookuptable->read_table("GSL_test_table.dat");
    std::cout << "Look up tau etc" << std::endl;
    start = clock();
    for (int i=0; i<100; i++){
//      lookuptable->compute(0, true);  // Look up tau etc
       lookuptable->compute(0, true, true);  // With root-polishing
    }
    duration = (double)(clock() - start) / CLOCKS_PER_SEC;
    std::cout << "Looked up tau  " << pmf->tau << " time " << duration/100. << std::endl;
    std::cout << pmf->tau << " " << pmf->dtaudf << " " << pmf->dtauds << std::endl;
        
    // Compute all possible conditional models
    int N = 10;
    double *eta = (double *) malloc(N * sizeof(double));
    double *p = (double *) malloc(N * sizeof(double));
    double *cv = (double *) malloc(3 * N * sizeof(double));
    for (int i=0; i<N; i++)
      eta[i] = (i + 1) * 1. / (N+1);    

    printvec(eta, N, "eta");
    pmf->counterflow(eta, N, p, N);
    printvec(p, N, " counterflow");
    pmf->CSD(eta, N, p, N);
    printvec(p, N, " CSD");
    pmf->PDF(eta, N, p, N);
    printvec(p, N, " PDF");
    pmf->Laplace(eta, N, p, N);
    printvec(p, N, " Laplace");
    pmf->CV(eta, N, cv, 3, N);
    printvec(cv, N, " CV[0]");
    std::cout << " counterflow " << pmf->counterflow(0.2) << std::endl; 
    std::cout << " CSD " << pmf->CSD(0.2) << std::endl;
    std::cout << " PDF " << pmf->PDF(0.2) << std::endl;
    std::cout << " Laplace " << pmf->Laplace(0.2) << std::endl;
    pmf->CV(0.2, cv, 3);
    std::cout << " CV[0] " << cv[0] << std::endl;
    start = clock();
    for (int i=0; i<100000; i++) pmf->CSD_verify(eta, N, p, N);
    duration = (double)(clock() - start) / CLOCKS_PER_SEC;
    std::cout << "Time CSD2 " << duration / 100000. << " value " << p[4] << endl;
    start = clock();
    for (int i=0; i<100000; i++) pmf->CSD(eta, N, p, N);
    duration = (double)(clock() - start) / CLOCKS_PER_SEC;
    std::cout << "Time CSD  " << duration / 100000. << " value " << p[4] << endl;
}
