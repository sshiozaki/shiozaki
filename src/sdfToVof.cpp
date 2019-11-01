/**
 * @file sdfToVof.cpp
 * @brief VOFC class
 * @author Tomohiro Otani
*/

#include "vof_conv.h"
using namespace std;

// #################################################################
/**
 * @brief ヘビサイド関数を用いたVOF作成
 * @param [in] Data SDF関数
 */
float* VOFC::SDFtoVOF(const float *Data){

  float *volumefraction = new float[subdinfo.dims[0]*subdinfo.dims[1]*subdinfo.dims[2]];

  unsigned long ret=0;
  float dx = (subdinfo.dx[0] + subdinfo.dx[1] + subdinfo.dx[2])/3e0;
  float tmp;

  for(unsigned int k=0;k<subdinfo.dims[2];k++){
    for(unsigned int j=0;j<subdinfo.dims[1];j++){
      for(unsigned int i=0;i<subdinfo.dims[0];i++){

        if(Data[ret]< -1e0*dx){
            volumefraction[ret] = 0e0;
        }else if(Data[ret]>dx){
            volumefraction[ret] = 1e0;
        }else{
            tmp = heavisideFunction(Data[ret],dx);
            if(tmp<0e0){
                volumefraction[ret] = 0e0;
            }else if(tmp>1e0){
                volumefraction[ret] = 1e0;
            }else{
                volumefraction[ret] = tmp;
            }
        }
          ret++;
      }
    }
  }
  return volumefraction;
}
// #################################################################
/**
 * @brief ヘビサイド関数 (fai+ep)/2pi + sin(pi*fai/ep) / 2pi
 * @param [in] fai f
 * @param [in] ep epsilon
 */
float VOFC::heavisideFunction(const float &fai,const float &ep){
  return 5e-1*(fai+ep)/(ep) + 5e-1*sin(M_PI*fai/ep)/M_PI;
}
