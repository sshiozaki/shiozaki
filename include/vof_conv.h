#ifndef _VOF_CONV_
#define _VOF_CONV_

//##################################################################################
//
// vof_conv.h
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
 * @file vof_conv.h
 * @brief VOFC class header
 * @author Tomohiro Otani
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <sys/stat.h>
#include "sdfInterface.h"
#include "info.h"
#include "vof_fileIO.h"

class VOFC{
 public:
  int numOfMesh;
  double *Data;
  float *vof, *sdf;
  meshInfo minfo;
  domainInfo dinfo;
  subDomainInfo subdinfo; //(MPI)
  VOFInfo vinfo;
  TextParser tp;

  std::string tmp_fname1, tmp_fname2;

 private:
  //SDFlib
  sdfPolyMesh mesh;
  sdfPtsWithNorm *pPwn;
  sdfSimpleVoxel *pSdf,*pSdf1,*pVof;

 public:
  VOFC(){};
  ~VOFC(){};

  void vof_main();
  void vof_main_2(const int count, const int countmax);

  void make_stlfile(const int count, const int countmax);
  void stl_read_binary(const std::string &fname, const int numoftriangle, float n1[][3], float x1[][3][3], 
                   unsigned char *buf, unsigned char *buf2);
  void stl_deformation(const int numoftriangle, const float n1[][3], const float x1[][3][3], float n2[][3], float x2[][3][3]);
  void stl_deformation_curve_3(const int numoftriangle, const float x1[][3][3], float n2[][3], float x2[][3][3], const int count, const int countmax);

  void crossProduct2(const double (&a)[3],const double (&b)[3], double c[3]);
  void stl_write_binary(const std::string &fname3, const int numoftriangle, const float n1[][3], const float x1[][3][3], 
                   const unsigned char *buf, const unsigned char *buf2);

  void initialize();
  void VOF_readGeometry();
  void VOF_readGeometry_2(const int count);

  void VOF_converter();
  void VOF_outputData(const int &numOfSubDomain);
  void VOF_outputData_2(const int &numOfSubDomain, const int count);
  void exportPVTI(const int num);

 private:
  //void make_stlfile(meshInfo &minfo, const int count);

  void read_mesh(sdfPolyMesh &mesh, const meshInfo &minfo);
  void read_mesh_2(sdfPolyMesh &mesh, meshInfo &minfo, const int count);

  float *SDFtoVOF(const float *Data);
  inline float heavisideFunction(const float &fai, const float &ep);
  float *copyArray(const float *DataArray, const size_t(&dim)[3],
                   const size_t(&dimStart)[3], const size_t(&dimEnd)[3]);
  // float *copyArray(const float *DataArray, const size_t(&dim)[3],const size_t(&dim2)[3]);
  double *copyArray(const float *DataArray, const size_t(&dim)[3],const size_t(&dim2)[3]);
  double *copyArray(float *DataArray, const size_t(&dim)[3]);
  void outputInfo2TextParser(const int &numOfSubDomain);
};

#endif //_VOF_CONV_
