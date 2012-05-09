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

#include "Integrator.h"

namespace pmfpack
{
  Integrator::Integrator(double **_parameters)
  : parameters(_parameters){}
  
  double function(const double phi, void *params)
  {
    double y, X, z;
    double **data = (double **) params;
    
    z = ((*data[3]) * phi - (*data[2])) / sqrt(1 - (*data[3]) * (*data[3]));
    X = erf(z);
    y = exp(-phi * phi / 2);
    return y * X * X / (*data[4]);
  }    
}
