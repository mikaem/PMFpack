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

#include "GSLRoot.h"

using namespace std;

namespace pmfpack{
    
  void my_error_handler (const char *reason, const char *file, int line, int gsl_errno)
  {
    gsl_stream_printf ("ERROR", file, line, reason);
    fflush (stdout);
    fprintf (stderr, "My GSL error handler invoked.\n");
    fflush (stderr);
    return;
  }

  GSLRoot::GSLRoot()
  : Root()
  {
    gsl_set_error_handler(&my_error_handler);
  }

  GSLRoot::GSLRoot(double _data[], Integrator *_integrator)
  : Root(_data, _integrator)
  {        
    gsl_set_error_handler(&my_error_handler);
  }
    
}
