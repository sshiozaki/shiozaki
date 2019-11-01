#ifndef _INFO_H_
#define _INFO_H_

//##################################################################################
//
// Info.h
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
 * @file info.h
 * @brief  domainInfo Class Header
 * @brief  meshInfo Class Header
 * @brief  VOFInfo Class Header
 * @author Tomohiro Otani
 */

#include <string>
#include "TextParser.h"
#include "sdfInterface.h"
#include "vof_define.h"

class domainInfo{
 public:
  domainInfo(){};
  ~domainInfo(){};

  int unitOfLength;
  float origin[3];    //origin of whole voxel
  float region[3];    //
  float dx[3];        //voxel size
  size_t numOfVoxel[3];    //number of voxel
  size_t MPIsubDivision[3]; //number of MPI subdomain
 void getDomainInfo(TextParser &tp);
 private:
};

class subDomainInfo{
  public:
  size_t MPIsubDivision[3];
  size_t dims[3];
  SDF::FVec3 origin;
  SDF::FVec3 dx;
  float *Data;
  void setSubDomain(const domainInfo &dinfo,const int &ix,const int &iy,const int &iz);
 private:
};

class meshInfo{
 public:
  meshInfo(){};
  ~meshInfo(){};
  std::string alias;
  std::string file, stlfile_name;
  int format;
  int type;
  int unitOfLength;
  void getMeshInfo(TextParser &tp);
 private:
};

class VOFInfo{
 public:
  VOFInfo(){};
  ~VOFInfo(){};
  bool reinitFMM;
  bool createVOF;
  bool outputSDF;
  bool outputVOF;
  bool insideOut;
  double isosurface;
  int type;
  int outputFormat;
  int visualizeFormat;
  std::string outputDir;
  std::string outputFilePath;
  std::string visualizeFilePath;
  void getVOFInfo(TextParser &tp);
 private:
};

#endif // _INFO_H__