// Copyright (C) 2012 Mikael Mortensen
//
// This file is part of PMFpack.
//
// PMFpack is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) a later version.
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
// First added:  2014-05-23
// Last changed: 2014-05-23

#include "HDF5Lookup.h"
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

namespace pmfpack
{
  
  HDF5Lookup::HDF5Lookup(Derivator *_derivator)
  : Lookup(_derivator) 
  {
    s = NULL;    
    fm = NULL;
    t = gsl_interp_polynomial;
    set_order(4);
  }  
  
  void HDF5Lookup::init(int N1, int N2)
  {
    Lookup::init(N1, N2);
    fm = (double *) malloc(Nf * sizeof(double));
    Is = (double *) malloc(Ns * sizeof(double));
    atangrid(1.0, 1.0);
    tau_table = f3tensor(5, Nf, Ns);
  } 
  
  HDF5Lookup::~HDF5Lookup()
  {
    lookup_clean();
  }
  
  void HDF5Lookup::lookup_clean()
  {
    free(fm);
    free(Is);
    free_f3tensor(tau_table);
    gsl_interp_free(s);
  }
  
  int HDF5Lookup::read_table(char* filename)
  {
    int status;
        
    H5File file = H5File(filename, H5F_ACC_RDONLY);
        
    DataSet fm_dat = file.openDataSet( "fmean" );
    DataSet is_dat = file.openDataSet( "Intensity of segregation" );
    hsize_t dims[1];
    fm_dat.getSpace().getSimpleExtentDims(dims);
    Nf = dims[0];
    is_dat.getSpace().getSimpleExtentDims(dims);
    Ns = dims[0];
    HDF5Lookup::init(Nf, Ns);
    fm_dat.read(fm, PredType::NATIVE_DOUBLE);
    is_dat.read(Is, PredType::NATIVE_DOUBLE);
    DataSet tau_dat = file.openDataSet( "Tau" );
    hsize_t Ndims[3] = {5, Nf, Ns};  // dataset dimensions
    DataSpace dataspace3D = tau_dat.getSpace();
    tau_dat.read(tau_table[0][0], PredType::NATIVE_DOUBLE, dataspace3D);
    Attribute attribute = tau_dat.openAttribute("Grid type");
    uint gt[1];
    attribute.read(PredType::NATIVE_INT, gt);
    grid_type = gt[0];
    std::cout << "Grid type " << grid_type << std::endl;

//     std::cout << " Read table" << std::endl;    
//     for (uint i=0; i<5; i++)
//     {
//       for (uint j=0; j<Nf; j++)
//       {
//         for (uint k=0; k<Ns; k++)
//           std::cout << i << " " << j << " " << k << " " << tau_table[i][j][k] << std::endl;
//       }
//     }

    
    return 1;
  }
  
  void HDF5Lookup::set_order(int _order)
  {
    if (s != NULL)
      gsl_interp_free(s);
    order = _order;
    s = gsl_interp_alloc(t, order);
  }
  
  void HDF5Lookup::generate_table(int N1, int N2, char *filename)
  {
    // Generate table using high precision
    
    if (fm != NULL && (N1 != Nf || N2 != Ns))
    {
        lookup_clean();
        init(N1, N2);
    }        
    else if (fm == NULL) 
    {
        init(N1, N2);
    }
    
    // Get high precision mean values
    typedef cpp_dec_float_24 T;

    T long_fm[Nf];
    T long_Is[Ns];
    const T pi = boost::math::constants::pi<T>();
    T at = atan2(pi, T(1.0));
    T dy = 2 * pi / (Ns + 1);
    for(int i = 0; i < Ns; i++)
      long_Is[i] = -pi + dy + dy * i;
    
    for(int i = 0; i < Ns; i++)
      long_Is[i] = (atan2(long_Is[i], T(1.0)) / at + 1) / 2;
    
    T dx = 2 * pi / (Nf + 1);
    for(int i = 0; i < Nf; i++)
      long_fm[i] = -pi + dx + dx * i;
    
    for(int i = 0; i < Nf; i++)
      long_fm[i] = (atan2(long_fm[i], T(1.0)) / at + 1) / 2;
    
    grid_type = 1;
    
    // Open a new file and dataset.
    H5File* file = new H5File(filename, H5F_ACC_TRUNC);
    
    // Stor fmean and Is
    hsize_t dim[1] = {Nf}; 
    DataSpace *dataspace = new DataSpace (1, dim);;
    DataSet *fm_dat = new DataSet(file->createDataSet( "fmean", 
                           PredType::NATIVE_DOUBLE, *dataspace) );
    
    dim[0] = Ns;
    dataspace->setExtentSimple(1, dim);
    DataSet *is_dat = new DataSet(file->createDataSet( "Intensity of segregation", 
                           PredType::NATIVE_DOUBLE, *dataspace) );
    
    fm_dat->write(fm, PredType::NATIVE_DOUBLE);
    is_dat->write(Is, PredType::NATIVE_DOUBLE);
    
    // Create a dataset for tau and its derivatives
    hsize_t dims[3] = {5, Nf, Ns};
    DataSpace *dataspace3D = new DataSpace (3, dims);
    DataSet *tau_dat = new DataSet(file->createDataSet( "Tau", 
                                   PredType::NATIVE_DOUBLE, *dataspace3D) );
    
    dim[0] = 1;
    uint gt[1] = {grid_type};
    DataSpace attr_dataspace = DataSpace (1, dim);
    Attribute attribute = tau_dat->createAttribute( "Grid type", PredType::NATIVE_UINT, 
                                                    attr_dataspace);
    attribute.write(PredType::NATIVE_UINT, gt);
    
    // Compute tau and its derivatives with high precision. 32 significant digits
    // are necessary to get 16 digits accuracy for second derivatives. 
    // Store using double - accurate to last digit.

    T h = T(1.0) / pow(10, 6);
    DerivativesHP<T>* der = new DerivativesHP<T>(T(0.5), T(0.1), h, 11, true, 2, 9);
    std::cout << std::scientific << std::setprecision(numeric_limits<T>::digits10);
    std::vector<T> result;
    for(int i = 0; i < Nf; i++){
      for(int j = 0; j < Ns; j++){
        T sigma = long_Is[j] * long_fm[i] * (1 - long_fm[i]);
        der->set_parameters(long_fm[i], sigma);
        result = der->derivatives();
        tau_table[0][i][j] = (double) result[0]; // tau       
        tau_table[1][i][j] = (double) result[1]; // dtaudf
        tau_table[2][i][j] = (double) result[3]; // dtauds
        tau_table[3][i][j] = (double) result[2]; // d2taudfdf
        tau_table[4][i][j] = (double) result[4]; // d2taudsds
        std::cout << i << " " << j << " "  << result[0] << " "<< long_fm[i] << " " << long_Is[j] << std::endl;
      }
    }
    tau_dat->write(tau_table[0][0], PredType::NATIVE_DOUBLE);
        
    file->close();
  }    
  
  void HDF5Lookup::gslpolin2(double x1a[], double x2a[], double **ya, double *y)
  {
    double *ymtmp = (double*) malloc(order * sizeof(double));
    for (int j = 0; j < order; j++) { 
      gsl_interp_init(s, x2a, ya[j], order);
      ymtmp[j] = gsl_interp_eval(s, x2a, ya[j], *sigma / *fmean / (1. - (*fmean)), NULL);
    }
    gsl_interp_init(s, x1a, ymtmp, order);
    *y = gsl_interp_eval(s, x1a, ymtmp, *fmean, NULL);
    free(ymtmp);
  }
  
  void HDF5Lookup::compute(int verbose, bool derivatives, bool polish)
  {
    uint m, n, ms1, ms2, ns1, ns2, ms, ns, check;
    double *x1s, *x2s, dy;

    // Locate point in lookuptable
    locate(&m, &n);
    
    // Compute if point is outside of table borders
    if (m == 0 || n == 0 || m == Nf || n == Ns)
    {
      if (derivatives)
        derivator->compute(verbose);
      else
        derivator->roots->froot->compute(verbose);
      return;
    }
    // Locate submatrix
    // Linear interpolation requires at least two points
    ms1 = IMIN(m+2, order), ms2 = IMAX(2, Nf-m+2), (ms1 < ms2)  ?  (ms=ms1) : (ms=ms2);
    ns1 = IMIN(n+2, order), ns2 = IMAX(2, Ns-n+2), (ns1 < ns2)  ?  (ns=ns1) : (ns=ns2);
    
    if (ms1 < order){
      ms2 = 2 * order - ms1;
    }
    else if(ms2 < 4){
      ms1 = 2 * order - ms2;
    }
    else{
      ms1 = order;
      ms2 = order;
    }
    if(ns1 < order){
      ns2 = 2 * order - ns1;
    }
    else if(ns2 < 4){
      ns1 = 2 * order - ns2;
    }
    else{
      ns1 = order;
      ns2 = order;
    }
    ms1 = (int)(floor(ms1 / 2.));
    ms2 = (int)(ceil(ms2 / 2.));
    ns1 = (int)(floor(ns1 / 2.));
    ns2 = (int)(ceil(ns2 / 2.));
    x1s = &fm[m-ms1];
    x2s = &Is[n-ns1];
    
    // Copy submatrix to ys and interpolate
    for(int i = m-ms1; i < m+ms2; i++){
      ys[i-m+ms1] = tau_table[0][i]+n-ns1;      
    }
    gslpolin2(x1s, x2s, ys, tau);
    // Polish the root with an fdf solver for best possible accuracy. The lookup is then initial guess.
    if (polish){
      if (derivatives)
        derivator->compute(verbose);
      else
        derivator->roots->fdfroot->compute(verbose);
      return; // Return because you don't need to look up derivatives you have already computed
    }
    if (derivatives){
      for (int k = 1; k < 5; k++){
        for(int i = m-ms1; i < m+ms2; i++){
          ys [i-m+ms1] = tau_table[k][i]+n-ns1;
        }
        gslpolin2(x1s, x2s, ys, dtau[k]);
      }    
    }
  }

  double ** HDF5Lookup::get_table(int i)
  {
    return tau_table[i];
  }
}
