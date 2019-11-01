/**
 * @file output.cpp
 * @brief VOFC classes
 * @author Tomohiro Otani
 */
#include "vof_conv.h"
#include <fstream>
#include <sys/stat.h>


using namespace std;

// #################################################################
/**
 * @brief データ出力
 */
void VOFC::VOF_outputData(const int &numOfSubDomain)
{
  string file,name1 = "VOF",name2 = "SDF";
  unsigned long ret;
  float origin[3],dx[3];
  size_t dims[3];

  for(int i=0;i<3;i++) dims[i]   = subdinfo.dims[i]+1;
  for(int i=0;i<3;i++) origin[i] = subdinfo.origin[i];
  for(int i=0;i<3;i++) dx[i]     = subdinfo.dx[i];

  mkdir(vinfo.outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  outputInfo2TextParser(numOfSubDomain);

  if(vinfo.outputVOF==ON){

    Data = new double[subdinfo.dims[0] * subdinfo.dims[1] * subdinfo.dims[2]];
    Data = copyArray(pVof->m_pData,dims,subdinfo.dims);

    for(int i=0;i<subdinfo.dims[0] * subdinfo.dims[1] * subdinfo.dims[2];i++){
      if(1.0e0<Data[i]) Data[i] = 1.0e0;
      if(0.0e0>Data[i]) Data[i] = 0.0e0;
    }

    if(vinfo.outputFormat==ORIGIN){
      file = vinfo.outputDir + "/chi_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + ".bin";
      vof_fileIO::exportVOF(file,Data,subdinfo.dims,vinfo.type);
    }
    if(vinfo.visualizeFormat==VTI){
      file = vinfo.outputDir + "/chi_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + ".vti";
      vof_fileIO::exportVTIBinary(file,name1,Data,subdinfo.dims,subdinfo.dims,dinfo.origin,dx,subdinfo.MPIsubDivision,volumeFraction);
    }
    delete[] Data;
  }

  if(vinfo.outputSDF==ON){
    Data = new double[dims[0] * dims[1] * dims[2]];

    if(vinfo.reinitFMM==ON)  Data = copyArray(pSdf1->m_pData,dims);
    if(vinfo.reinitFMM==OFF) Data = copyArray(pSdf->m_pData, dims);

    // if(vinfo.outputFormat==HDF5){
    //   file = vinfo.outputDir + "/SDF_" + to_string(subdinfo.MPIsubDivision[0]) + "_"
    //                                    + to_string(subdinfo.MPIsubDivision[1]) + "_"
    //                                    + to_string(subdinfo.MPIsubDivision[2]) + "_" + vinfo.outputFilePath;
    //   fileIO::export_HDF5(file,name2,Data,dims,origin,dx);
    // }
      file = vinfo.outputDir + "/sdf_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + ".vti";
      // cout << file << endl;
      vof_fileIO::exportVTIBinary(file,name2,Data,subdinfo.dims,dims,dinfo.origin,dx,subdinfo.MPIsubDivision,signedDistanceFunction);
       file = vinfo.outputDir + "/sdf_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                        + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                        + to_string(subdinfo.MPIsubDivision[2]) + ".bin";
      vof_fileIO::exportVOF(file,Data,dims,vinfo.type);
    delete[] Data;
  }
}

void VOFC::VOF_outputData_2(const int &numOfSubDomain, const int count)
{
  string file,name1 = "VOF",name2 = "SDF";
  unsigned long ret;
  float origin[3],dx[3];
  size_t dims[3];

  for(int i=0;i<3;i++) dims[i]   = subdinfo.dims[i]+1;
  for(int i=0;i<3;i++) origin[i] = subdinfo.origin[i];
  for(int i=0;i<3;i++) dx[i]     = subdinfo.dx[i];

  mkdir(vinfo.outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  outputInfo2TextParser(numOfSubDomain);

  if(vinfo.outputVOF==ON){

    Data = new double[subdinfo.dims[0] * subdinfo.dims[1] * subdinfo.dims[2]];
    Data = copyArray(pVof->m_pData,dims,subdinfo.dims);

    for(int i=0;i<subdinfo.dims[0] * subdinfo.dims[1] * subdinfo.dims[2];i++){
      if(1.0e0<Data[i]) Data[i] = 1.0e0;
      if(0.0e0>Data[i]) Data[i] = 0.0e0;
    }

    if(vinfo.outputFormat==ORIGIN){
      file = vinfo.outputDir + "/chi_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + ".bin";
      vof_fileIO::exportVOF(file,Data,subdinfo.dims,vinfo.type);
    }
    if(vinfo.visualizeFormat==VTI){
      file = vinfo.outputDir + "/chi_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + ".vti";
      vof_fileIO::exportVTIBinary(file,name1,Data,subdinfo.dims,subdinfo.dims,dinfo.origin,dx,subdinfo.MPIsubDivision,volumeFraction);
    }
    delete[] Data;
  }

  if(vinfo.outputSDF==ON){
    Data = new double[dims[0] * dims[1] * dims[2]];

    if(vinfo.reinitFMM==ON)  Data = copyArray(pSdf1->m_pData,dims);
    if(vinfo.reinitFMM==OFF) Data = copyArray(pSdf->m_pData, dims);

    // if(vinfo.outputFormat==HDF5){
    //   file = vinfo.outputDir + "/SDF_" + to_string(subdinfo.MPIsubDivision[0]) + "_"
    //                                    + to_string(subdinfo.MPIsubDivision[1]) + "_"
    //                                    + to_string(subdinfo.MPIsubDivision[2]) + "_" + vinfo.outputFilePath;
    //   fileIO::export_HDF5(file,name2,Data,dims,origin,dx);
    // }
      file = vinfo.outputDir + "/sdf_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                       + to_string(subdinfo.MPIsubDivision[2]) + "_" + to_string(count) + ".vti";
      // cout << file << endl;
      vof_fileIO::exportVTIBinary(file,name2,Data,subdinfo.dims,dims,dinfo.origin,dx,subdinfo.MPIsubDivision,signedDistanceFunction);
       file = vinfo.outputDir + "/sdf_" + to_string(subdinfo.MPIsubDivision[0]) + "-"
                                        + to_string(subdinfo.MPIsubDivision[1]) + "-"
                                        + to_string(subdinfo.MPIsubDivision[2]) + ".bin";
      vof_fileIO::exportVOF(file,Data,dims,vinfo.type);
    delete[] Data;
  }
}

// #################################################################
/**
 * @brief tpへのregionデータ書き込み(MPI対応)
 * @param [in] numOfSubDomain MPI domain name
 */
void VOFC::outputInfo2TextParser(const int &numOfSubDomain)
{
  string label, str;

  label = "/region" + to_string(numOfSubDomain) + "/NumOfVoxel";
  str = "( " + to_string(subdinfo.dims[0]) + ", " + to_string(subdinfo.dims[1]) + ", " + to_string(subdinfo.dims[2]) + ")";
  tp.createLeaf(label, str);

  label = "/region" + to_string(numOfSubDomain) + "/GlobalLocation";
  str = "( " + to_string(subdinfo.MPIsubDivision[0]) + ", "
             + to_string(subdinfo.MPIsubDivision[1]) + ", "
             + to_string(subdinfo.MPIsubDivision[2]) + ")";
  tp.createLeaf(label,str);

  label = "/region" + to_string(numOfSubDomain) + "/GlobalVoxelOrigin";
  str = "( " + to_string(subdinfo.origin[0]) + ", "
             + to_string(subdinfo.origin[1]) + ", "
             + to_string(subdinfo.origin[2]) + ")";
  tp.createLeaf(label,str);
}


// #################################################################
/**
 * @brief 配列の部分コピー
 * @param [in] DataArray 配列
 * @param [in] dim voxel数
 */
double *VOFC::copyArray(float *DataArray, const size_t(&dim)[3])
{
  double *copyDataArray = new double[(dim[0] * dim[1] * dim[2])];
  unsigned long ret = 0;

  for (unsigned int k=0;k<dim[2];k++) {
    for (unsigned int j=0;j<dim[1];j++){
      for (unsigned int i=0;i<dim[0];i++){
        copyDataArray[ret] = DataArray[i + j*dim[0] + k*dim[0]*dim[1]];
        ret++;
       }
    }
  }
  return copyDataArray;
}

// #################################################################
/**
 * @brief 配列の部分コピー
 * @param [in] DataArray 配列
 * @param [in] dim voxel数
 * @param [in] dimStart voxel始点
 * @param [in] dimEnd voxel終点
 */
double *VOFC::copyArray(const float *DataArray, const size_t(&dim)[3],const size_t(&dim2)[3])
{
  double *copyDataArray = new double[dim2[0]*dim2[1]*dim2[2]];
  unsigned long ret = 0;

  for (auto k=0;k<dim2[2];k++) {
    for (auto j=0;j<dim2[1];j++){
      for (auto i=0;i<dim2[0];i++){
        copyDataArray[ret] = (double)DataArray[i + j*dim[0] + k*dim[0]*dim[1]];
        ret++;
       }
    }
  }
  return copyDataArray;
}

// #################################################################
/**
 * @brief 配列の部分コピー
 * @param [in] DataArray 配列
 * @param [in] dim voxel数
 * @param [in] dimStart voxel始点
 * @param [in] dimEnd voxel終点
 */
float *VOFC::copyArray(const float *DataArray, const size_t(&dim)[3],
                       const size_t(&dimStart)[3], const size_t(&dimEnd)[3]) {
  float *copyDataArray =
      new float[(dimEnd[0] - dimStart[0]) * (dimEnd[1] - dimStart[1]) * (dimEnd[2] - dimStart[2])];
  unsigned long ret = 0;

  for (auto k=dimStart[2];k<dimEnd[2];k++) {
    for (auto j=dimStart[1];j<dimEnd[1];j++){
      for (auto i=dimStart[0];i<dimEnd[0];i++){
        copyDataArray[ret] = DataArray[i + j*dim[0] + k*dim[0]*dim[1]];
        ret++;
       }
    }
  }
  return copyDataArray;
}
// #################################################################
/**
 * @brief export parallel VTI file
 * @param [in] dataname
 */
void VOFC::exportPVTI(const int num)
{
  fstream ofs;
  string dataName;

  if(num==volumeFraction) dataName = "chi";
  if(num==signedDistanceFunction) dataName = "sdf";
  string file = vinfo.outputDir + "/" + "total_" + dataName + ".pvti";

  cout << file << endl;

  ofs.open(file.c_str(),ios::out);

  ofs << "<?xml version=\"1.0\"?>" << endl;

  ofs << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
  ofs << "  <PImageData WholeExtent=\" 0 " << dinfo.numOfVoxel[0] << " 0 " << dinfo.numOfVoxel[1] << " 0 " << dinfo.numOfVoxel[2]
         << "\" GhostLevel=\"0\" Origin=\"" << dinfo.origin[0] << " " << dinfo.origin[1] << " " << dinfo.origin[2]
         << "\" Spacing=\"" << dinfo.dx[0] << " " << dinfo.dx[1] << " " << dinfo.dx[2] << "\">" << endl;

  if(num==signedDistanceFunction){
    ofs << "<PPointData Scalars=\""<< dataName << "\">" << endl;
    ofs << "<PDataArray type=\"Float64\" Name=\"" << dataName <<"\"/>" << endl;
    ofs << "</PPointData>" << endl;
    ofs << "<PCellData>" << endl;
    ofs << "</PCellData>" << endl;
  }else if(num==volumeFraction){
    ofs << "<PPointData>" << endl;
    ofs << "</PPointData>" << endl;
    ofs << "<PCellData Scalars=\""<< dataName << "\">" << endl;
    ofs << "<PDataArray type=\"Float64\" Name=\"" << dataName <<"\"/>" << endl;
    ofs << "</PCellData>" << endl;
  }
  for(unsigned int iz=0;iz<dinfo.MPIsubDivision[2];iz++){
    for(unsigned int iy=0;iy<dinfo.MPIsubDivision[1];iy++){
      for(unsigned int ix=0;ix<dinfo.MPIsubDivision[0];ix++){
        ofs << "<Piece Extent=\"" << subdinfo.dims[0]*ix << " " << subdinfo.dims[0]*(ix+1) << " "
                                  << subdinfo.dims[1]*iy << " " << subdinfo.dims[1]*(iy+1) << " "
                                  << subdinfo.dims[2]*iz << " " << subdinfo.dims[2]*(iz+1)
             << "\" Source=\""<< dataName <<"_"<< ix << "-" << iy << "-" << iz << ".vti" << "\"/>" << endl;
      }
    }
  }
  ofs<<"</PImageData>"<<endl;
  ofs<<"</VTKFile>"<<endl;

}
