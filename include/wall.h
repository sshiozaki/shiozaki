//##################################################################################
//
// Coil deployment simulator (Otani et al., Med. & Biol. Eng. & Compt., 2017)
//
// Copyright (c) 2012 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################
#ifndef _WALL_H_
#define _WALL_H_

/**
 * @file   wall.h
 * @brief  coil deployment main Header
 * @author T. Otani
 */
#include "coRotationalBeam_define.h"
#include "TextParser.h"

class Wall{
 public:
  double origin[3];
  double glboalLength[3];
  int numOfVoxel[3];
  double dx[3];
  DOUBLEARRAY3 SDF;
  std::string SDFfile;

  void inputParameters(TextParser &tp);
  void inputParameters_catheter(TextParser &tp);
  void C3D8_N(double (&N)[8],const double &g1,const double &g2,const double &g3);
  void C3D8_dNdr(double (&dNdr)[8][3],const double &g1,const double &g2,const double &g3);
  void exportVTIBinary(const std::string &file,const std::string &name);
 private:
  void readSDF(const std::string &file,const int &nx,const int &ny,const int &nz);
};

#endif