#ifndef _VOF_FILEIO_H_
#define _VOF_FILEIO_H_

//##################################################################################
//
// fileIO.h
//
// Copyright (c) 2016 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file fileIO.h
 * @brief  fileIO Class Header
 * @author Tomohiro Otani
 */

#include <string>
#include "vof_define.h"

class vof_fileIO{
 public:

  static void exportVTKAscii(const std::string &file,const std::string &name,const float *Data,
                            const size_t (&dim)[3],const float (&ori)[3],const float (&aspectRatio)[3]);
  static void exportVTIAscii(const std::string &file,const std::string &name,const float *Data,
                            const size_t (&dim)[3],const float (&ori)[3],const float (&dx)[3]);

  static void exportVTIBinary(const std::string &file,const std::string &name,const double *Data,
                              const size_t (&dim)[3],const float (&ori)[3],const float (&dx)[3]);

  static void exportVTIBinary(const std::string &file,const std::string &name,const double *Data,const size_t (&dim)[3],
                             const size_t (&dim_data)[3],const float (&ori)[3],const float (&dx)[3],const size_t (&MPIsubDivision)[3],const int num);

  static void exportVOF(const std::string &file,const double *Data,const size_t (&size)[3],const int &type);

};

#endif // _VOF_FILEIO_H__