//##################################################################################
//
// Multiple beam element simulator
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################

#ifndef _FILEIO_H_
#define _FILEIO_H_

/**
 * @file   fileIO.h
 * @brief  original file input Header
 * @author T. Otani
 */

#include <string>
#include "coRotationalBeam_define.h"

class fileIO{
 public:
 static int CountNumbersOfTextLines(const std::string &filePath );
 static void read_original_file_double(DOUBLEARRAY2 &x,const int &nd,const int &number,const std::string &file);
 static void read_original_file_int(INTARRAY2 &x,const int &nd,const int &number,const std::string &file);
};
#endif //_IO_ORIGINAL_FILE_HPP_