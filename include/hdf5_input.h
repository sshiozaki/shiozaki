#ifndef _HDF5_INPUT_H_
#define _HDF5_INPUT_H_

//##################################################################################
//
// hdf5_input.h
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   hdf5_input.h
 * @brief  hdf5 input Header
 * @author T.Otani
 */
#include <string>
#include <H5Cpp.h>
#include "allocation.h"

class HDF5input{
 public:
  static void read_1dINTARRAY(H5::H5File &file,const std::string &dataName,INTARRAY1 &data);
  static void read_2dINTARRAY(H5::H5File &file,const std::string &dataName,INTARRAY2 &data);
  static void read_3dINTARRAY(H5::H5File &file,const std::string &dataName,INTARRAY3 &data);
  static void read_4dINTARRAY(H5::H5File &file,const std::string &dataName,INTARRAY4 &data);

  static void read_1dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,DOUBLEARRAY1 &data);
  static void read_2dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,DOUBLEARRAY2 &data);
  static void read_3dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,DOUBLEARRAY3 &data);
  static void read_4dDOUBLEARRAY(H5::H5File &file,const std::string &dataName,DOUBLEARRAY4 &data);
};

#endif
