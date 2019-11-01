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

/**
 * @file wall.cpp
 * @author T. Otani
 */
#include "wall.h"
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


// #################################################################
/**
 * @brief shape function of voxel
 * @param [out] N   shape function
 * @param [in] g1   position in x direction
 * @param [in] g2   position in y direction
 * @param [in] g3   position in z direction
 */
void Wall::C3D8_N(double (&N)[8],const double &g1,const double &g2,const double &g3)
{
  N[0] = 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0-g3);
  N[1] = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0-g3);
  N[2] = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0-g3);
  N[3] = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0-g3);
  N[4] = 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0+g3);
  N[5] = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0+g3);
  N[6] = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0+g3);
  N[7] = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0+g3);
}

// #################################################################
/**
 * @brief derivative of shape function of voxel
 * @param [out] dNdr   derivative of shape function
 * @param [in] g1      position in x direction
 * @param [in] g2      position in y direction
 * @param [in] g3      position in z direction
 */
void Wall::C3D8_dNdr(double (&dNdr)[8][3],const double &g1,const double &g2,const double &g3)
{
  dNdr[0][0] = -1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr[0][1] = -1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr[0][2] = -1.25e-1 * (1e0-g1) * (1e0-g2);

  dNdr[1][0] =  1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr[1][1] = -1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr[1][2] = -1.25e-1 * (1e0+g1) * (1e0-g2);

  dNdr[2][0] =  1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr[2][1] =  1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr[2][2] = -1.25e-1 * (1e0+g1) * (1e0+g2);

  dNdr[3][0] = -1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr[3][1] =  1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr[3][2] = -1.25e-1 * (1e0-g1) * (1e0+g2);

  dNdr[4][0] = -1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr[4][1] = -1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr[4][2] =  1.25e-1 * (1e0-g1) * (1e0-g2);

  dNdr[5][0] =  1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr[5][1] = -1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr[5][2] =  1.25e-1 * (1e0+g1) * (1e0-g2);

  dNdr[6][0] =  1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr[6][1] =  1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr[6][2] =  1.25e-1 * (1e0+g1) * (1e0+g2);

  dNdr[7][0] = -1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr[7][1] =  1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr[7][2] =  1.25e-1 * (1e0-g1) * (1e0+g2);
}

// #################################################################
/**
 * @brief read SDF from file (binary format)
 * @param [in] file     SDF file name
 * @param [in] nx       number of voxel (x)
 * @param [in] ny       number of voxel (y)
 * @param [in] nz       number of voxel (z)
 */
void Wall::readSDF(const string &file,const int &nx,const int &ny,const int &nz)
{
  ifstream fin;
  double *data  = (double *)malloc(nx*ny*nz*sizeof(double));

  fin.open(file.c_str(),ios::in|ios::binary);
  if(!fin){
    cout << "invalid file name-> " << file << endl;
    exit(0);
  }

  unsigned long allsize = sizeof(double)*nx*ny*nz;

  fin.read((char *)data,allsize);
  fin.close();

  // #pragma omp parallel for
  for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
       SDF[i][j][k] = data[i + j*nx + k*nx*ny];
     }
    }
  }

  free(data);
}

// #########################################################
/**
 * @brief input parameters from tp file
 * @param [in]    tp   TextParser file
 */
void Wall::inputParameters(TextParser &tp)
{
  string str,base_label,label;

  // Global_origin
  label = "/wallDomain/GlobalOrigin";
  if ( !tp.getInspectedVector(label, origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_region
  label = "/wallDomain/GlobalLength";
  if ( !tp.getInspectedVector(label, glboalLength, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_voxel
  label = "/wallDomain/GlobalVoxel";
  if ( !tp.getInspectedVector(label, numOfVoxel, 3) ){
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // sdf dir
  label = "/wallDomain/sdfFile";
  if ( !tp.getInspectedValue(label, SDFfile) ){
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  for(int i=0;i<3;i++) dx[i] = glboalLength[i]/(double)numOfVoxel[i];
  SDF=Allocation::allocate3dDOUBLE(numOfVoxel[0]+1,numOfVoxel[1]+1,numOfVoxel[2]+1);
  readSDF(SDFfile,numOfVoxel[0]+1,numOfVoxel[1]+1,numOfVoxel[2]+1);
  string file="test.vti";
  string name="SDF";
  exportVTIBinary(file,name);
}

// #########################################################
/**
 * @brief input parameters from tp file
 * @param [in]    tp   TextParser file
 */
void Wall::inputParameters_catheter(TextParser &tp)
{
  string str,base_label,label;

  // Global_origin
  label = "/wallDomain/GlobalOrigin";
  if ( !tp.getInspectedVector(label, origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_region
  label = "/wallDomain/GlobalLength";
  if ( !tp.getInspectedVector(label, glboalLength, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // Global_voxel
  label = "/wallDomain/GlobalVoxel";
  if ( !tp.getInspectedVector(label, numOfVoxel, 3) ){
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  // sdf dir
  label = "/wallDomain_catheter/sdfFile";
  if ( !tp.getInspectedValue(label, SDFfile) ){
    cout << "ERROR : in parsing [" << label << "]" << endl;
    exit(0);
  }

  for(int i=0;i<3;i++) dx[i] = glboalLength[i]/(double)numOfVoxel[i];
  SDF=Allocation::allocate3dDOUBLE(numOfVoxel[0]+1,numOfVoxel[1]+1,numOfVoxel[2]+1);
  readSDF(SDFfile,numOfVoxel[0]+1,numOfVoxel[1]+1,numOfVoxel[2]+1);
  string file="test_catheter.vti";
  string name="SDF_catheter";
  exportVTIBinary(file,name);
}

// #################################################################
/**
 * @brief vti(binary) export
 * @param [in] file vtk file name
 * @param [in] name output data name
 */
void Wall::exportVTIBinary(const string &file,const string &name)
{
  fstream ofs;

  ofs.open(file.c_str(),ios::out);
  //header書き込み
  ofs << "<?xml version=\"1.0\"?>"<<endl;
  ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
  ofs << "  <ImageData WholeExtent= \"0 " << numOfVoxel[0] << " 0 " << numOfVoxel[1] << " 0 " << numOfVoxel[2] <<
                          "\" Origin= \""   << origin[0] << " "   << origin[1] << " "   << origin[2] <<
                          "\" Spacing= \""   << dx[0]  << " "   << dx[1]  << " "   << dx[2]  << "\">" << endl;
  ofs << "    <Piece Extent= \"0 " << numOfVoxel[0] << " 0 " << numOfVoxel[1] << " 0 " << numOfVoxel[2] << "\" > " << endl;
  ofs << "      <PointData Scalars=\""<< name.c_str() << "\">" << endl;
  ofs << "        <DataArray type=\"Float64\" Name=\""<< name.c_str() << "\" format=\"appended\" offset=\"0\">" << endl;
  ofs << "        </DataArray>" << endl;
  ofs << "      </PointData>" <<endl;
  ofs << "    </Piece>" <<endl;
  ofs << "  </ImageData>" << endl;
  ofs << "  <AppendedData encoding=\"raw\">" << endl;
  ofs << "    _";
  ofs.close();

  //データ書き込み
  ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
  unsigned long allsize=(numOfVoxel[0]+1)*(numOfVoxel[1]+1)*(numOfVoxel[2]+1)*sizeof(double);
  double *data=(double *)malloc(allsize*sizeof(double));
  int num=0;
  for(int k=0;k<numOfVoxel[2]+1;k++){
    for(int j=0;j<numOfVoxel[1]+1;j++){
      for(int i=0;i<numOfVoxel[0]+1;i++){
        data[num]=SDF[i][j][k];
        num++;
      }
    }
  }

  ofs.write((char *)&allsize,sizeof(allsize));
  ofs.write((char*)data,allsize);
  ofs.close();

  ofs.open(file.c_str(),ios::out | ios::app);
  ofs << endl;
  ofs << "  </AppendedData>" << endl;
  ofs << "</VTKFile>" << endl;
  free(data);
}