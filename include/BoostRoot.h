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

#ifndef __BOOSTROOT_H
#define __BOOSTROOT_H

#include <boost/math/tools/roots.hpp>
#include <boost/function.hpp>
#include "Root.h"

namespace pmfpack
{
  
  class BoostRoot : public Root
  {
  public:
    
    BoostRoot(Integrator *);
    
    virtual ~BoostRoot() {};
            
    virtual double compute(int) = 0;
    
    virtual void realloc(const unsigned int) = 0;
    
    virtual void error_message_off() {};
        
  };    
}

#endif
