
/**
 * @file fileIO.cpp
 * @brief  domainInfo Class Header
 * @brief  meshInfo Class Header
 * @brief  VOFInfo Class Header
 * @author Tomohiro Otani
 */

#include "vof_fileIO.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

// #################################################################
/**
 * @brief vtk(ascii)出力 (debug用、データの大規模化に伴い使用を控える)
 * @param [in] file vtiファイル名
 * @param [in] name 出力データ名
 * @param [in] dim ボクセル数(X,Y,Z)
 * @param [in] Data data格納ポインタ
 * @param [in] ori 原点座標(X,Y,Z)
 * @param [in] aspectRatio ボクセル比率(X,Y,Z)
 */
void vof_fileIO::exportVTKAscii(const string &file,const string &name,const float *Data,
                            const size_t (&dim)[3],const float (&ori)[3],const float (&aspectRatio)[3])
{
  ofstream ofs(file.c_str());

  ofs << "# vtk DataFile Version 2.0" << endl;
  ofs << "test" << endl;
  ofs << "ASCII" << endl;
  ofs << "DATASET STRUCTURED_POINTS" << endl;
  ofs << "DIMENSIONS " << dim[0]+1<<" "<< dim[1]+1 << " " << dim[2]+1 << endl;
  ofs << "ORIGIN " << ori[0]<<" "<< ori[1] << " " << ori[2] << endl;
  ofs << "ASPECT_RATIO " << aspectRatio[0]<<" " << aspectRatio[1] << " " << aspectRatio[2] << endl;
  ofs << endl;
  ofs << "CELL_DATA " <<dim[0]*dim[1]*dim[2] << endl;
    ofs << "SCALARS "<< name.c_str() <<" float " << endl;
    ofs << "LOOKUP_TABLE default" << endl;

  int ret=0;
  for(unsigned int k=0;k<dim[2];k++){
    for(unsigned int j=0;j<dim[1];j++){
      for(unsigned int i=0;i<dim[0];i++){
          ofs << Data[ret] <<" ";
          ret++;
          ofs << endl;
      }
    }
  }

}
// #################################################################
/**
 * @brief vti(ascii)出力 (debug用、データの大規模化に伴い使用を控える)
 * @param [in] file vtiファイル名
 * @param [in] name 出力データ名
 * @param [in] dim ボクセル数(X,Y,Z)
 * @param [in] Data data格納ポインタ
 * @param [in] ori 原点座標(X,Y,Z)
 * @param [in] dx ボクセルサイズ(X,Y,Z)
 */
void vof_fileIO::exportVTIAscii(const string &file,const string &name,const float *Data,
                            const size_t (&dim)[3],const float (&ori)[3],const float (&dx)[3])
{
  ofstream ofs(file.c_str());

  ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
  ofs << "  <ImageData WholeExtent= \"0 " << dim[0] << " 0 " << dim[1] << " 0 " << dim[2] <<
                          "\" Origin= \""   << ori[0] << " "   << ori[1] << " "   << ori[2] <<
                          "\" Spacing= \""   << dx[0]  << " "   << dx[1]  << " "   << dx[2]  << "\">" << endl;
  ofs << "    <Piece Extent= \"0 " << dim[0] << " 0 " << dim[1] << " 0 " << dim[2] << "\" > " << endl;
  ofs << "      <PointData>" << endl;
  ofs << "      </PointData>" << endl;
  ofs << "      <CellData Scalars=\""<< name.c_str() << "\">" << endl;
  ofs << "        <DataArray type=\"Float32\" Name=\""<< name.c_str() <<"\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1.0\">" << endl;

  int ret=0;
  for(unsigned int k=0;k<dim[2];k++){
    for(unsigned int j=0;j<dim[1];j++){
      for(unsigned int i=0;i<dim[0];i++){
          ofs << "        " << Data[ret] <<" ";
          ret++;
          if(ret%6==0) ofs << endl;
      }
    }
  }

  ofs << "      </DataArray>" << endl;
  ofs << "    </CellData>" <<endl;
  ofs << "  </Piece>" <<endl;
  ofs << "  </ImageData>" << endl;
  ofs << "</VTKFile>" << endl;
}

// #################################################################
/**
 * @brief vti(binary)出力
 * @param [in] file vtkファイル名
 * @param [in] name 出力データ名
 * @param [in] dim ボクセル数(X,Y,Z)
 * @param [in] Data data格納ポインタ
 * @param [in] ori 原点座標(X,Y,Z)
 * @param [in] dx ボクセルサイズ(X,Y,Z)
 */
void vof_fileIO::exportVTIBinary(const string &file,const string &name,const double *Data,
                             const size_t (&dim)[3],const float (&ori)[3],const float (&dx)[3])
{
  fstream ofs;

  ofs.open(file.c_str(),ios::out);

  //header
  ofs << "<?xml version=\"1.0\"?>"<<endl;
  ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
  ofs << "  <ImageData WholeExtent= \"" << dim[0] << " " << dim[0] << " "
                                        << dim[1] << " " << dim[1] << " "
                                        << dim[2] << " " << dim[2] <<
                          "\" Origin= \""   << ori[0] << " "   << ori[1] << " "   << ori[2] <<
                          "\" Spacing= \""   << dx[0]  << " "   << dx[1]  << " "   << dx[2]  << "\">" << endl;
  ofs << "    <Piece Extent= \""  << dim[0] << " " << dim[0] << " "
                                  << dim[1] << " " << dim[1] << " "
                                  << dim[2] << " " << dim[2] << "\">" << endl;

  // ofs << "      <PointData>" << endl;
  // ofs << "      </PointData>" << endl;
  // ofs << "      <CellData Scalars=\""<< name.c_str() << "\">" << endl;
  // ofs << "        <DataArray type=\"Float32\" Name=\""<< name.c_str()
  //           << "\" format=\"appended\" RangeMin=\""<< 0 <<"\" RangeMax=\""<< 1 <<"\" offset=\"0\">" << endl;
  // ofs << "        </DataArray>" << endl;
  // ofs << "      </CellData>" <<endl;

  ofs << "      <PointData Scalars=\""<< name.c_str() << "\">" << endl;
  ofs << "        <DataArray type=\"Float64\" Name=\""<< name.c_str()
            << "\" format=\"appended\" RangeMin=\""<< 0 <<"\" RangeMax=\""<< 1 <<"\" offset=\"0\">" << endl;
  ofs << "        </DataArray>" << endl;
  ofs << "      </PointData>" << endl;
  ofs << "      <CellData>" << endl;
  ofs << "      </CellData>" <<endl;

  ofs << "    </Piece>" <<endl;
  ofs << "  </ImageData>" << endl;
  ofs << "  <AppendedData encoding=\"raw\">" << endl;
  ofs << "    _";
  ofs.close();

  //data
  ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
  unsigned long allsize=(dim[0]+1)*(dim[1]+1)*(dim[2]+1)*sizeof(double);
  ofs.write((char *)&allsize,sizeof(allsize));
  ofs.write((char *)Data,allsize);
  ofs.close();

  ofs.open(file.c_str(),ios::out | ios::app);
  ofs << endl;
  ofs << "  </AppendedData>" << endl;
  ofs << "</VTKFile>" << endl;
}
// #################################################################
/**
 * @brief vti(binary)出力
 * @param [in] file vtkファイル名
 * @param [in] name 出力データ名
 * @param [in] dim ボクセル数(X,Y,Z)
 * @param [in] Data data格納ポインタ
 * @param [in] ori 原点座標(X,Y,Z)
 * @param [in] dx ボクセルサイズ(X,Y,Z)
 * @param [in] MPIsubDivision subdomain(X,Y,Z)
 */
void vof_fileIO::exportVTIBinary(const string &file,const string &name,const double *Data,const size_t (&dim)[3],
                             const size_t (&dim_data)[3],const float (&ori)[3],const float (&dx)[3],const size_t (&MPIsubDivision)[3],const int num)
{
  fstream ofs;

  ofs.open(file.c_str(),ios::out);

  //header
  ofs << "<?xml version=\"1.0\"?>"<<endl;
  ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
  ofs << "  <ImageData WholeExtent= \"" << MPIsubDivision[0] * dim[0] << " " << (MPIsubDivision[0]+1) * dim[0] << " "
                                        << MPIsubDivision[1] * dim[1] << " " << (MPIsubDivision[1]+1) * dim[1] << " "
                                        << MPIsubDivision[2] * dim[2] << " " << (MPIsubDivision[2]+1) * dim[2] <<
                          "\" Origin= \""   << ori[0] << " "   << ori[1] << " "   << ori[2] <<
                          "\" Spacing= \""   << dx[0]  << " "   << dx[1]  << " "   << dx[2]  << "\">" << endl;
  ofs << "    <Piece Extent= \""  << MPIsubDivision[0] * dim[0] << " " << (MPIsubDivision[0]+1) * dim[0] << " "
                                  << MPIsubDivision[1] * dim[1] << " " << (MPIsubDivision[1]+1) * dim[1] << " "
                                  << MPIsubDivision[2] * dim[2] << " " << (MPIsubDivision[2]+1) * dim[2] << "\">" << endl;

  if(num==signedDistanceFunction){
  ofs << "      <PointData Scalars=\""<< name.c_str() << "\">" << endl;
  ofs << "        <DataArray type=\"Float64\" Name=\""<< name.c_str()
            << "\" format=\"appended\" RangeMin=\""<< 0 <<"\" RangeMax=\""<< 1 <<"\" offset=\"0\">" << endl;
  ofs << "        </DataArray>" << endl;
  ofs << "      </PointData>" << endl;
  ofs << "      <CellData>" << endl;
  ofs << "      </CellData>" <<endl;
 }else if(num==volumeFraction){
  ofs << "      <PointData>" << endl;
  ofs << "      </PointData>" << endl;
  ofs << "      <CellData Scalars=\""<< name.c_str() << "\">" << endl;
  ofs << "        <DataArray type=\"Float64\" Name=\""<< name.c_str()
            << "\" format=\"appended\" RangeMin=\""<< 0 <<"\" RangeMax=\""<< 1 <<"\" offset=\"0\">" << endl;
  ofs << "        </DataArray>" << endl;
  ofs << "      </CellData>" <<endl;
 }

  ofs << "    </Piece>" <<endl;
  ofs << "  </ImageData>" << endl;
  ofs << "  <AppendedData encoding=\"raw\">" << endl;
  ofs << "    _";
  ofs.close();

  //data
  ofs.open(file.c_str(),ios::out | ios::app | ios_base::binary);
  unsigned long allsize=dim_data[0]*dim_data[1]*dim_data[2]*sizeof(double);
  ofs.write((char *)&allsize,sizeof(allsize));
  ofs.write((char *)Data,allsize);
  ofs.close();


  ofs.open(file.c_str(),ios::out | ios::app);
  ofs << endl;
  ofs << "  </AppendedData>" << endl;
  ofs << "</VTKFile>" << endl;
}

// #################################################################
/**
 * @brief vof(binary)出力
 * @param [in] file vtkファイル名
 * @param [in] name 出力データ名
 * @param [in] oinfo output設定構造体
 */
void vof_fileIO::exportVOF(const string &file,const double *Data,const size_t (&size)[3],const int &type)
{
  fstream ofs;

  //データ書き込み
  if(type==BINARY){
    ofs.open(file.c_str(),ios::out | ios_base::binary);
    unsigned long allsize=size[0]*size[1]*size[2]*sizeof(double);
    ofs.write((char*)Data,allsize);
    ofs.close();
  }

  if(type==ASCII){
    int num=0;
    ofs.open(file.c_str(),ios::out);
    for(int iz=0;iz<size[2];iz++){
      for(int iy=0;iy<size[1];iy++){
        for(int ix=0;ix<size[0];ix++){
          ofs << ix << " " << iy << " " << iz << " " << Data[num] << endl;
          num++;
        }
      }
    }
  }

}


