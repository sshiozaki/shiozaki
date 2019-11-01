/**
 * @file info.cpp
 * @brief meshinfo class
 * @brief domaininfo class
 * @brief VOFinfo class
 * @author Tomohiro Otani
 */
#include "info.h"
#include <iostream>
#include <cstdio>

using namespace std;

// #################################################################
/**
 * @brief VOF設定情報読み込み
 * @param [in] tp TextParser用ファイル
 */
void VOFInfo::getVOFInfo(TextParser &tp){

  string str;
  string label;

  label = "/VOF/reinitFMM";
  if ( !tp.getInspectedValue(label,str) ) exit(0);
  if(!strcasecmp(str.c_str(),"ON"))  reinitFMM = ON;
  if(!strcasecmp(str.c_str(),"OFF")) reinitFMM = OFF;

  label = "/VOF/createVOF";
  if ( !tp.getInspectedValue(label,str) ) exit(0);
  if(!strcasecmp(str.c_str(),"ON"))  createVOF = ON;
  if(!strcasecmp(str.c_str(),"OFF")) createVOF = OFF;

  label = "/VOF/outputSDF";
  if ( !tp.getInspectedValue(label,str) ) exit(0);
  if(!strcasecmp(str.c_str(),"ON"))  outputSDF = ON;
  if(!strcasecmp(str.c_str(),"OFF")) outputSDF = OFF;

  label = "/VOF/outputVOF";
  if ( !tp.getInspectedValue(label,str) ) exit(0);
  if(!strcasecmp(str.c_str(),"ON"))  outputVOF = ON;
  if(!strcasecmp(str.c_str(),"OFF")) outputVOF = OFF;

  if(outputSDF==OFF && outputVOF==OFF){
    cout << "no output file. please confirm\n";
    exit(0);
  }

  label = "/VOF/outputDir";
  if ( !tp.getInspectedValue(label,outputDir)){
    cout << "Output directory is not set" << endl;
    exit(0);
  }

  label = "/VOF/type";
  if ( !tp.getInspectedValue(label,str)){
    cout << "Output directory is not set" << endl;
    exit(0);
  }
  if(!strcasecmp(str.c_str(),"ASCII"))  type = ASCII;
  if(!strcasecmp(str.c_str(),"BINARY")) type = BINARY;

  label = "/VOF/outputFormat";
  if ( !tp.getInspectedValue(label,str)){
    cout << label << "not found" << endl;
    exit(1);
  }
  if(!strcasecmp(str.c_str(),"hdf5"))  outputFormat = HDF5;
  if(!strcasecmp(str.c_str(),"VOF"))   outputFormat = ORIGIN;

  label = "/VOF/visualizeFormat";
  if ( !tp.getInspectedValue(label,str)){
    cout << label << " not found" << endl;
    exit(1);
  }
  if(!strcasecmp(str.c_str(),"vti"))  visualizeFormat = VTI;

  label = "/VOF/isosurface";
  if ( !tp.getInspectedValue(label,isosurface)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  label = "/VOF/insideOut";
  if ( !tp.getInspectedValue(label,str) ){
    cout << label << " not found" << endl;
    exit(1);
  }
  if(!strcasecmp(str.c_str(),"ON"))  insideOut = ON;
  if(!strcasecmp(str.c_str(),"OFF")) insideOut = OFF;
}

// #################################################################
/**
 * @brief ドメイン情報読み込み
 * @param [in] tp TextParserファイル
 */
void domainInfo::getDomainInfo(TextParser &tp){

  int tmp[3];
  string str;
  string label,label_m,label_leaf;

  // Global_origin
  label = "/DomainInfo/GlobalOrigin";
  if ( !tp.getInspectedVector(label, origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_region
  label = "/DomainInfo/GlobalRegion";
  if ( !tp.getInspectedVector(label, region, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_voxel
  label = "/DomainInfo/GlobalVoxel";
  if ( !tp.getInspectedVector(label, tmp, 3) ){
    cout << label << " is not set" << endl;
    exit(0);
  }
  for(int i=0;i<3;i++) numOfVoxel[i] = tmp[i];

   // Global_division
  label = "/DomainInfo/GlobalDivision";
  if ( !tp.getInspectedVector(label, tmp, 3) ){
    cout << label << " is not set" << endl;
    exit(0);
  }
  for(int i=0;i<3;i++) MPIsubDivision[i] = tmp[i];

  for(int i=0;i<3;i++) dx[i] = region[i]/(double)numOfVoxel[i];
}

// #################################################################
/**
 * @brief ドメイン情報読み込み
 * @param [in] tp TextParserファイル
 */
void subDomainInfo::setSubDomain(const domainInfo &dinfo,const int &ix,const int &iy,const int &iz)
{
  for(int i=0;i<3;i++){
    origin[i] = dinfo.origin[i];
    dx[i] = dinfo.dx[i];
    dims[i] = dinfo.numOfVoxel[i] / dinfo.MPIsubDivision[i];
  }

  MPIsubDivision[0] = ix;
  MPIsubDivision[1] = iy;
  MPIsubDivision[2] = iz;

  origin[0] += dx[0] * ix * dims[0];
  origin[1] += dx[1] * iy * dims[1];
  origin[2] += dx[2] * iz * dims[2];

}
// #################################################################
/**
 * @brief mesh情報読み込み
 * @param [in] mesh mesh情報設定用クラス
 * @param [in] m meshの数
 */
void meshInfo::getMeshInfo(TextParser &tp){

  string str,label;

  label = "/Mesh/filePath";
  if ( !tp.getInspectedValue(label,file)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  label = "/Mesh/format";
  if ( !tp.getInspectedValue(label,str)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  if(!strcasecmp(str.c_str(),"STL")) format = STL;
  if(!strcasecmp(str.c_str(),"OBJ")) format = OBJ;


  label = "/Mesh/type";
  if ( !tp.getInspectedValue(label,str)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  if(!strcasecmp(str.c_str(),"ASCII"))  type = ASCII;
  if(!strcasecmp(str.c_str(),"BINARY")) type = BINARY;

  label = "/Mesh/filename";
  if ( !tp.getInspectedValue(label,stlfile_name)){
    cout << label << " is not set" << endl;
    exit(0);
  }
}
