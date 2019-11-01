/**
 * @file vof_converter.cpp
 * @brief VOFC class
 * @author Tomohiro Otani
 */
#include "vof_conv.h"

using namespace std;

// #################################################################
/**
 * @brief SDF, VOF作成, 単一meshファイル用
 */
void VOFC::VOF_readGeometry()
{
  double tm1,tm2;

  //load mesh
    tm1 = SDF::GetTime();
    read_mesh(mesh,minfo);
    tm2 = SDF::GetTime();
    printf("load geometry: %ld points, %ld faces | %lf sec.\n",mesh.getNumPts(), mesh.getNumFcs(), tm2-tm1);

  //get points+normal data from mesh
  tm1 = SDF::GetTime();
  pPwn = SDF::GetPWNfromMESH(mesh);
  if(!pPwn){
    fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
    exit(0);
  }
  tm2 = SDF::GetTime();
  printf("Convert %ld points of pointSet | %f sec. \n",pPwn->getNumPts(),tm2-tm1);

}


// #################################################################
/**
 * @brief SDF, VOF作成, 単一meshファイル用
 */
void VOFC::VOF_readGeometry_2(const int count)
{
  double tm1,tm2;

  //load mesh
    tm1 = SDF::GetTime();
    read_mesh_2(mesh,minfo, count);
    tm2 = SDF::GetTime();
    printf("load geometry: %ld points, %ld faces | %lf sec.\n",mesh.getNumPts(), mesh.getNumFcs(), tm2-tm1);

  //get points+normal data from mesh
  tm1 = SDF::GetTime();
  pPwn = SDF::GetPWNfromMESH(mesh);
  if(!pPwn){
    fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
    exit(0);
  }
  tm2 = SDF::GetTime();
  printf("Convert %ld points of pointSet | %f sec. \n",pPwn->getNumPts(),tm2-tm1);

}
// #################################################################
/**
 * @brief for single mesh only
 */
void VOFC::VOF_converter(){

  double tm1,tm2;

  for(int i=0;i<3;i++) subdinfo.dims[i] += 1;

  //get SDF(simple voxel) from points+normals
  tm1 = SDF::GetTime();
  pSdf = SDF::GetSDFfromPWN(*pPwn,subdinfo.dims,subdinfo.origin,subdinfo.dx);
  if(!pSdf){
    fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
    exit(0);
  }
  tm2 = SDF::GetTime();
  printf("create Sv SDF: dims(%ld,%ld,%ld) | %f sec.\n",subdinfo.dims[0],subdinfo.dims[1],subdinfo.dims[2],tm2-tm1);

  //reinit SDF value by using FMM scheme
  if(vinfo.reinitFMM==ON){
    tm1 = SDF::GetTime();
    pSdf1 = SDF::ReinitSDF_FMM(*pSdf);
    if(!pSdf){
      fprintf(stderr, "%s\n",SDF::s_lastErrMsg.c_str());
      exit(0);
    }
    tm2 = SDF::GetTime();
    printf("reinit Sv SDF (FMM scheme) | %f sec. \n",tm2-tm1);
  }

  int tmp=0;
  for(unsigned int k=0;k<subdinfo.dims[2];k++){
    for(unsigned int j=0;j<subdinfo.dims[1];j++){
      for(unsigned int i=0;i<subdinfo.dims[0];i++){
        if(vinfo.reinitFMM==ON)  pSdf1->m_pData[tmp] += vinfo.isosurface;
        if(vinfo.reinitFMM==OFF) pSdf->m_pData[tmp]  += vinfo.isosurface;
        tmp++;
      }
    }
  }

 //rasterlize: get VOF from SDF
  if(vinfo.createVOF==ON){
    tm1 = SDF::GetTime();
    vof = new float[subdinfo.dims[0]*subdinfo.dims[1]*subdinfo.dims[2]];
      //if(vinfo.reinitFMM==ON)  vof = SDFtoVOF(pSdf1->m_pData);
      //if(vinfo.reinitFMM==OFF) vof = SDFtoVOF(pSdf->m_pData);
      if(vinfo.reinitFMM==ON)  pVof = SDF::GetVOFfromSDF(*pSdf1);
      if(vinfo.reinitFMM==OFF) pVof = SDF::GetVOFfromSDF(*pSdf);
    if(!vof){
      fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
      exit(0);
    }
    tm2 = SDF::GetTime();
    printf("create Sv VOF(Heaviside function): dims(%ld, %ld, %ld) | %f sec.\n",subdinfo.dims[0],subdinfo.dims[1],subdinfo.dims[2],tm2-tm1);
  }

  if(vinfo.insideOut==ON){
    int tmp=0;
    for(unsigned int k=0;k<subdinfo.dims[2];k++){
      for(unsigned int j=0;j<subdinfo.dims[1];j++){
        for(unsigned int i=0;i<subdinfo.dims[0];i++){
          pVof->m_pData[tmp] = 1.0e0-pVof->m_pData[tmp];
          tmp++;
        }
      }
    }
  }

  for(int i=0;i<3;i++) subdinfo.dims[i] -= 1;

}
// #########################################################
/**
 * @brief read mesh file
 * @param [in,out] mesh SDFlib mesh class
 * @param [in] minfo mesh info. class
 */
void VOFC::read_mesh(sdfPolyMesh &mesh,const meshInfo &minfo){

  switch (minfo.format){
    case STL:
      if(minfo.type==ASCII){
        if(!mesh.loadSla(minfo.file)){
          fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
          exit(0);
        }
      }
      if(minfo.type==BINARY){
        if(!mesh.loadSlb(minfo.file)){
          fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
          exit(0);
        }
      }
      break;
    case OBJ:
      if(!mesh.loadObj(minfo.file)){
        fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
        exit(0);
      }
      break;
    default:
      cout << "Invalid file format" <<endl;
      exit(1);
  }
}

// #########################################################
/**
 * @brief read mesh file
 * @param [in,out] mesh SDFlib mesh class
 * @param [in] minfo mesh info. class
 */
void VOFC::read_mesh_2(sdfPolyMesh &mesh,meshInfo &minfo, const int count){
  // if(count > 1){
  //   minfo.file = minfo.stlfile_name + "_" + to_string(count) + ".stl";
  // }

  minfo.file = minfo.stlfile_name + "_" + to_string(count) + ".stl";
  printf("%d\n", count);

  switch (minfo.format){
    case STL:
      if(minfo.type==ASCII){
        if(!mesh.loadSla(minfo.file)){
          fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
          exit(0);
        }
      }
      if(minfo.type==BINARY){
        if(!mesh.loadSlb(minfo.file)){
          fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
          exit(0);
        }
      }
      break;
    case OBJ:
      if(!mesh.loadObj(minfo.file)){
        fprintf(stderr,"%s\n",SDF::s_lastErrMsg.c_str());
        exit(0);
      }
      break;
    default:
      cout << "Invalid file format" <<endl;
      exit(1);
  }
}

// ################################################
/**
 * @brief initilaize and input various info.
 */
void VOFC::initialize(){

  //meshファイル情報読み込み
  minfo.getMeshInfo(tp);

  //ドメイン情報読み込み
//  domainInfo dinfo;
  dinfo.getDomainInfo(tp);

  //VOF作成情報読み込み
//  VOFInfo vinfo;
  vinfo.getVOFInfo(tp);
}


// ################################################
/**
 * @brief initilaize and input various info.
 */
void VOFC::make_stlfile(const int count, const int countmax)
{
  unsigned char buf[80];
  unsigned char buf2[2];
  unsigned int numoftriangle;

  std::string read_stlfile_name, write_stlfile_name;

  if(count == 1){
    read_stlfile_name = minfo.file;
    write_stlfile_name = minfo.stlfile_name + "_" + to_string(count) + ".stl";
  }else{
    read_stlfile_name = minfo.stlfile_name + "_" + to_string(count-1) + ".stl";
    write_stlfile_name = minfo.stlfile_name + "_" + to_string(count) + ".stl";
  }


  FILE *fp;
  fp=fopen(read_stlfile_name.c_str(), "rb");
  fread(buf, 80, 1, fp);
  fread(&numoftriangle, 4, 1, fp);
  fclose(fp);

  float n1[numoftriangle][3];
  float x1[numoftriangle][3][3];
  float n2[numoftriangle][3];
  float x2[numoftriangle][3][3];

  stl_read_binary(read_stlfile_name, numoftriangle, n1, x1, buf, buf2);

  //stl_deformation(numoftriangle, n1, x1, n2, x2);
  stl_deformation_curve_3(numoftriangle, x1, n2, x2, count, countmax);
  stl_write_binary(write_stlfile_name, numoftriangle, n2, x2, buf, buf2);
  // if(count%5 == 0){

  //   std::string delete_file;
  //   for(int i=1; i<5; i++){
  //     delete_file = minfo.stlfile_name + "_" + to_string(count-i) + ".stl";
  //     remove(delete_file.c_str());
  //   }
  // }
}


void VOFC::crossProduct2(const double (&a)[3],const double (&b)[3], double c[3])
{
  double dc = 0e0;

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  for(int j=0;j<3;j++) dc += (c[j]*c[j]);

  dc = sqrt(dc);

  for(int j=0; j<3; j++){
    c[j] = c[j]/dc;
  }
}

void VOFC::stl_read_binary(const std::string &fname, const int numoftriangle, float n1[][3], float x1[][3][3], 
                   unsigned char *buf, unsigned char *buf2)
{
  int tmp;

  FILE *fp;
  fp=fopen(fname.c_str(), "rb");
  fread(buf, 80, 1, fp);
  fread(&tmp, 4, 1, fp);

  for(int i=0; i<numoftriangle; i++){
    for(int j=0; j<3; j++){
      fread(&n1[i][j], 4, 1, fp);
    }
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
        fread(&x1[i][j][k], 4, 1, fp);
      }
    }
    fread(buf2, 2, 1, fp);
  }
  fclose(fp);
}

void VOFC::stl_write_binary(const std::string &fname3, const int numoftriangle, const float n1[][3], const float x1[][3][3], 
                   const unsigned char *buf, const unsigned char *buf2)
{
  FILE *fp3;
  fp3 = fopen(fname3.c_str(), "wb");
  fwrite(&buf, 80, 1, fp3);
  fwrite(&numoftriangle, 4, 1, fp3);

  for(int i=0; i<numoftriangle; i++){
    for(int j=0; j<3; j++){
      fwrite(&n1[i][j], 4, 1, fp3);
    }
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
        fwrite(&x1[i][j][k], 4, 1, fp3);
      }
    }
    fwrite(buf2, 2, 1, fp3);
  }
  fclose(fp3);
}

void VOFC::stl_deformation(const int numoftriangle, const float n1[][3], const float x1[][3][3], float n2[][3], float x2[][3][3])
{
  double r, r2;
  double theta;
  double tmp1[3], tmp2[3], tmp3[3];

  for(int i=0; i<numoftriangle; i++){
    for(int j=0; j<3; j++){
      r = sqrt(x1[i][j][0]*x1[i][j][0]+x1[i][j][1]*x1[i][j][1]);
      theta = atan(x1[i][j][1]/x1[i][j][0]);
      if(x1[i][j][0]>0e0){
        theta = atan(x1[i][j][1]/x1[i][j][0]);
      }else{
        theta = atan(x1[i][j][1]/x1[i][j][0]);
        theta = theta + acos(-1);
      }

      r2 = r - 5e-2;

      x2[i][j][0] = r2*cos(theta);
      x2[i][j][1] = r2*sin(theta);
      x2[i][j][2] = x1[i][j][2];       
    }

    for(int k=0; k<3; k++){
      tmp1[k] = x2[i][1][k] - x2[i][0][k];
      tmp2[k] = x2[i][2][k] - x2[i][1][k];
    }
    crossProduct2(tmp1,tmp2,tmp3);

    for(int j=0; j<3; j++){
      n2[i][j] = tmp3[j];
    }
  }
}


void VOFC::stl_deformation_curve_3(const int numoftriangle, const float x1[][3][3], float n2[][3], float x2[][3][3], const int count, const int countmax)
{
  double r, r2;
  double theta;
  double tmp1[3], tmp2[3], tmp3[3];
  double tmp4, tmp5;
  double R_vessel=4e0;

  double R_vessel_2;

  if(count == 1){
      R_vessel_2 = 1000e0;
    }else if(count == 2){
      R_vessel_2 = 500e0;
    }else if(count == 3){
      R_vessel_2 = 400e0;
    }else if(count == 4){
      R_vessel_2 = 300e0;
    }else if(count == 5){
      R_vessel_2 = 200e0;
    }else if(count == 6){
      R_vessel_2 = 151e0 + R_vessel - 2e0;
    }else if(count == 7){
      R_vessel_2 = 131e0 + R_vessel - 2e0;
    }else if(count == 8){
      R_vessel_2 = 111e0 + R_vessel - 2e0;
    }else if(count == 9){
      R_vessel_2 = 102e0 + R_vessel - 2e0;
    }else if(count > 9 && count < 51){
    R_vessel_2 = 103e0 - sqrt(99e0*99e0 - ((double)count-9e0 - 100e0)*((double)count-9e0 - 100e0)); 
  }else if(count > 50 && count < 111){
    R_vessel_2 = 103e0 - sqrt(99e0*99e0 - (40e0 + ((double)count-9e0-40e0)*5e-1 - 100e0)*(40e0 + ((double)count-9e0-40e0)*5e-1 - 100e0));
  }else{
    R_vessel_2 = 103e0 - sqrt(99e0*99e0 - (70e0 + ((double)count-9e0-100e0)*3e-1 - 100e0)*(70e0 + ((double)count-9e0-100e0)*3e-1 - 100e0));
  }

  for(int i=0; i<numoftriangle; i++){
    for(int j=0; j<3; j++){
      tmp5 = x1[i][j][2] - 10e0;
      if(tmp5 > 0e0){
        if(tmp5 < (R_vessel*acos(-1e0)*5e-1)){
          for(int n1=0; n1<15710; n1++){
            theta = 1e-4 * (double) n1;
            if((R_vessel_2*theta) > tmp5){
              break;
            }
          }
          r = R_vessel_2 - x1[i][j][0];

          x2[i][j][0] = R_vessel_2 - r*cos(theta);
          x2[i][j][1] = x1[i][j][1];
          x2[i][j][2] = r*sin(theta) + 10e0;
        }else{
          for(int n1=0; n1<15710; n1++){
            theta = 1e-4 * (double) n1;
            if((R_vessel_2*theta) > (R_vessel*acos(-1e0)*5e-1)){
              break;
            }
          }
          r = R_vessel_2 - x1[i][j][0];
          tmp4 = tmp5 - (R_vessel*acos(-1e0)*5e-1);

          x2[i][j][0] = R_vessel_2 - r*cos(theta) + tmp4*sin(theta);
          x2[i][j][1] = x1[i][j][1];
          x2[i][j][2] = r*sin(theta) + 10e0 + tmp4*cos(theta);
        }
      }else{
        tmp5 = fabs(tmp5);
        if(tmp5 < (R_vessel*acos(-1e0)*5e-1)){
          for(int n1=0; n1<15710; n1++){
            theta = 1e-4 * (double) n1;
            if((R_vessel_2*theta) > tmp5){
              break;
            }
          }
          r = R_vessel_2 - x1[i][j][0];

          x2[i][j][0] = R_vessel_2 - r*cos(theta);
          x2[i][j][1] = x1[i][j][1];
          x2[i][j][2] = (-1e0)*r*sin(theta) + 10e0;
        }else{
          for(int n1=0; n1<157100; n1++){
            theta = 1e-5 * (double) n1;
            if((R_vessel_2*theta) > (R_vessel*acos(-1e0)*5e-1)){
              break;
            }
          }
          r = R_vessel_2 - x1[i][j][0];
          tmp4 = tmp5 - (R_vessel*acos(-1e0)*5e-1);

          x2[i][j][0] = R_vessel_2 + - r*cos(theta) + tmp4*sin(theta);
          x2[i][j][1] = x1[i][j][1];
          x2[i][j][2] = (-1e0)*(r*sin(theta) + tmp4*cos(theta)) + 10e0;
        }
      }
    }

    for(int k=0; k<3; k++){
      tmp1[k] = x2[i][1][k] - x2[i][0][k];
      tmp2[k] = x2[i][2][k] - x2[i][1][k];
    }
    crossProduct2(tmp1,tmp2,tmp3);

    for(int j=0; j<3; j++){
      n2[i][j] = tmp3[j];
    }   
  }
}
