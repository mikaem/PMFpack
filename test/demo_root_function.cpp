
#include "PMF.h"

using namespace pmfpack;
using namespace std;

int main()
{
    double fm;
    fm = 0.70837758;
    PMF *pmf = new PMF(0.70837758, 0.194360748 * fm * (1 - fm));
    
    pmf->set_fmean_gradient(-0.3496497, 0, 0);
    pmf->set_sigma_gradient(-0.4294126, 0, 0);
    pmf->chi = 1.003770;
    pmf->DT = 1.;

//    pmf->root->realloc(1);
//    pmf->initial_guess(0.2);
//    pmf->froot->error_message_off();
    pmf->verbose = 1;
    pmf->compute(0, true);
      
}