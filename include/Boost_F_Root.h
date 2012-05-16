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

#ifndef __BOOST_F_ROOT_H
#define __BOOST_F_ROOT_H

#include "BoostRoot.h"

namespace pmfpack
{
  struct BoostFunction
  {
    BoostFunction(void *a, double (*f)(double, void*)) 
    : params(a), func(f) {};
    double operator ()(double const& z) 
      {
        return func(z, params);        
      };
  private:
    boost::function<double (double, void*)> func;
    void *params;
  };
  
  class tolerance {
  public:
    tolerance(double eps) :
      _eps(eps) {
    }
    bool operator()(double a, double b) {
      return (fabs(b - a) <= _eps);
    }
  private:
    double _eps;
  };

  class Boost_F_Root : public BoostRoot
  {
  public:
    Boost_F_Root();
    
    Boost_F_Root(Integrator *);
    
    ~Boost_F_Root();
            
    BoostFunction *F;
    
    tolerance *tol;
    
    int fsolver;
    
    boost::uintmax_t maxiter;
    
    std::pair<double, double> found;
        
    double compute(int);
    
    void realloc(const unsigned int);
        
    double upper, lower;
    
  };  
}

#endif
