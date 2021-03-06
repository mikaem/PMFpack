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

%module PMFpack
%{
#define SWIG_FILE_WITH_INIT
#include "numpy/arrayobject.h"
#include "PMF.h"
using namespace pmfpack;
%}

%include "numpy.i"

%init %{
  import_array();
%}

%feature("autodoc",1);

%include "carrays.i"
%array_functions(double,double);
%array_functions(int,int);
%array_functions(double*,doublep);

%apply (double* INPLACE_ARRAY1, int DIM1) {(double*, int ), (double*, int)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double*, int, int)};
%include PMF.h
%include Lookup.h
%include GSLLookup.h
%include PMFMath.h

