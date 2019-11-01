#ifndef _HDF5_OUTPUT_H_
#define _HDF5_OUTPUT_H_

//##################################################################################
//
// hdf5_output.h
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   hdf5_output.h
 * @brief  hdf5 output Header
 * @author T.Otani
 */
#include <string>
#include <H5Cpp.h>
#include "allocation.h"

class HDF5output{
 public:
  static void write_1dINTARRAY(H5::H5File &file,const std::string &dataName,const INTARRAY1 &data,const int size0);
  static void write_2dINTARRAY(H5::H5File &file,const std::string &dataName,const INTARRAY2 &data,const int size0,const int size1);
  static void write_3dINTARRAY(H5::H5File &file,const std::string &dataName,const INTARRAY3 &data,const int size0,const int size1,const int size2);
  static void write_4dINTARRAY(H5::H5File &file,const std::string &dataName,const INTARRAY4 &data,const int size0,const int size1,const int size2,const int size3);

  static void write_1dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,const DOUBLEARRAY1 &data,const int size0);
  static void write_2dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,const DOUBLEARRAY2 &data,const int size0,const int size1);
  static void write_3dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,const DOUBLEARRAY3 &data,const int size0,const int size1,const int size2);
  static void write_4dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,const DOUBLEARRAY4 &data,const int size0,const int size1,const int size2,const int size3);
};

#endif
