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

#ifndef __BOOST_FDF_ROOT_H
#define __BOOST_FDF_ROOT_H

#include "BoostRoot.h"

namespace pmfpack
{
  struct BoostF2Function
  {
    BoostF2Function(Integrator *a, void (*f)(double, Integrator*, double*, double*)) 
    : integrator(a), func(f) {};
    boost::math::tuple<double, double> operator ()(double const& z) 
      {
        double ff, dff;
        func(z, integrator, &ff, &dff);
        return boost::math::make_tuple(ff, dff);
      };
  private:
    boost::function<void (double, Integrator*, double*, double*)> func;
    Integrator *integrator;
  };

  struct BoostF3Function
  {
    BoostF3Function(Integrator *a, void (*f)(double, Integrator*, double*, double*, double*)) 
    : integrator(a), func(f) {};
    boost::math::tuple<double, double, double> operator ()(double const& z) 
      {
        double ff, dff, df2;
        func(z, integrator, &ff, &dff, &df2);
        return boost::math::make_tuple(ff, dff, df2);
      };
  private:
    boost::function<void (double, Integrator*, double*, double*, double*)> func;
    Integrator *integrator;
  };
  
  class Boost_FDF_Root : public BoostRoot
  {
  public:
    
    Boost_FDF_Root(Integrator *);
    
    ~Boost_FDF_Root();
            
    BoostF2Function *F2;
    
    BoostF3Function *F3;
    
    int fdfsolver;
    
    int digits;
    
    void set_digits(const uint);
    
    boost::uintmax_t maxiter;
    
    double compute(int);
    
    void realloc(const unsigned int);
        
    double upper, lower;
    
  };  
}

#endif
