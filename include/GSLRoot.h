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

#ifndef __GSLROOT_H
#define __GSLROOT_H

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "Root.h"

namespace pmfpack
{
  void my_error_handler (const char *, const char *, int, int);
  double gslfunc(double x, void *params);
  double gsldfunc(double x, void *params);
  void gslfdfunc(double x, void *params, double *f, double *df);

  class GSLRoot : public Root
  {
  public:
    
    GSLRoot(Integrator *);
    
    virtual ~GSLRoot() {};
            
    virtual double compute(int) = 0;
    
    virtual void realloc(const unsigned int) = 0;
    
    virtual void error_message_off() {gsl_set_error_handler_off();};
    
  };    
}

#endif
