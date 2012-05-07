
#include "PMF.h"

using namespace pmfpack;
using namespace std;

int main()
{
    PMF *pmf = new PMF(0.5, 0.1);
    
    cout << pmf->data[0] << endl;
    
    cout << Erf(0.5) << " " << gsl_cdf_ugaussian_P(0.5) << endl; 
    cout << Erfinv(0.4) << " " << gsl_cdf_ugaussian_Pinv(0.4) << endl; 
//    pmf->root->realloc(1);
//    pmf->initial_guess(0.2);
    pmf->froot->error_message_off();
    pmf->verbose = 1;
    pmf->compute_tau(0);
//     pmf->compute_tau_and_derivatives();
    
        
}