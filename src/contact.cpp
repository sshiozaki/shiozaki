//##################################################################################
//
// beam deployment simulator
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################

/**
 * @file contact.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

// #########################################################
/**
 * @brief calc self contact force
 * @param [in]     ibeam    target beam number
 */
void multipleBeamSimulator::selfContact_APC(const int ibeam)
{
  int ia[3],ic2;
  double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
             +(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0));
  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;

  //#pragma omp parallel for private(ia,ic2,vcon)
  for(int ic1=0;ic1<beam[ibeam].numOfElm;ic1++){

    for(int i=0;i<3;i++) ia[i]=box[ibeam][ic1][i];
    if(ia[0]<0 || ia[0]>=numOfVoxel[0]) continue;
    if(ia[1]<0 || ia[1]>=numOfVoxel[1]) continue;
    if(ia[2]<0 || ia[2]>=numOfVoxel[2]) continue;

    beamContactTMP[ibeam][0][ic1].component.clear();

    for(int ix=ia[0]-1;ix<=ia[0]+1;ix++){
      if(ix<0 || ix>=numOfVoxel[0]) continue;
      for(int iy=ia[1]-1;iy<=ia[1]+1;iy++){
        if(iy<0 || iy>=numOfVoxel[1]) continue;
        for(int iz=ia[2]-1;iz<=ia[2]+1;iz++){
          if(iz<0 || iz>=numOfVoxel[2]) continue;

          for(int it=0;it<numOfNodeInBox[ix][iy][iz][ibeam];it++){

            ic2 = nodeInBox[ix][iy][iz][ibeam][it];
            if(ic2-ic1<2) continue;   //to avoid self contact and duplicated count
            ContactInBox(beamContactTMP[ibeam][0][ic1].component,ic1,ic2,beam[ibeam],beam[ibeam],beamContactHistory[ibeam][0][ic1].contactHistoryMAP,kcon,E_astalisk);
          }
        }
      }
    }
  }

  //assembly
  int icc;
  for(int ic1=0;ic1<beam[ibeam].numOfElm;ic1++){
    for(int ic2=0;ic2<beamContactTMP[ibeam][0][ic1].component.size();ic2++){
      icc=beamContactTMP[ibeam][ibeam][ic1].component[ic2].ic;
      for(int j=0;j<3;j++){
        beam[ibeam].fcon[beam[ibeam].ie[ic1][0]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].fn[0][j];
        beam[ibeam].fcon[beam[ibeam].ie[ic1][1]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].fn[1][j];
        beam[ibeam].fcon[beam[ibeam].ie[icc][0]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].fn[2][j];
        beam[ibeam].fcon[beam[ibeam].ie[icc][1]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].fn[3][j];
        beam[ibeam].fcon[beam[ibeam].ie[ic1][0]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].ft[0][j];
        beam[ibeam].fcon[beam[ibeam].ie[ic1][1]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].ft[1][j];
        beam[ibeam].fcon[beam[ibeam].ie[icc][0]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].ft[2][j];
        beam[ibeam].fcon[beam[ibeam].ie[icc][1]][j]+=beamContactTMP[ibeam][0][ic1].component[ic2].ft[3][j];
      }
    }
  }

}

// #########################################################
/**
 * @brief contact force [ibeam1 x ibeam2]
 * @param [in]     ibeam1    target beam number
 * @param [in]     ibeam2    target beam number (ibeam2>ibeam1)
 */
void multipleBeamSimulator::contact_beams_APC(const int ibeam1,const int ibeam2)
{
  int ia[3],ic2;
  double Poisson=5e-1*(beam[ibeam1].Poisson+beam[ibeam2].Poisson);
  double EA=5e-1*(beam[ibeam1].EA+beam[ibeam2].EA);
  double rad=5e-1*(beam[ibeam1].rad+beam[ibeam2].rad);
  double kcon=(1e0-Poisson*Poisson)/(EA*PI*rad*rad)+(1e0-Poisson*Poisson)/(EA*PI*rad*rad);
  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;

  //#pragma omp parallel for private(ia,ic2,vcon)
  for(int ic1=0;ic1<beam[ibeam1].numOfElm;ic1++){
  //for(int ic1=1;ic1<beam[ibeam1].numOfElm-1;ic1++){    //except 0 and N-1 node !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for(int i=0;i<3;i++) ia[i]=box[ibeam1][ic1][i];
    if(ia[0]<0 || ia[0]>=numOfVoxel[0]) continue;
    if(ia[1]<0 || ia[1]>=numOfVoxel[1]) continue;
    if(ia[2]<0 || ia[2]>=numOfVoxel[2]) continue;

    beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component.clear();

    for(int ix=ia[0]-1;ix<=ia[0]+1;ix++){
      if(ix<0 || ix>=numOfVoxel[0]) continue;
      for(int iy=ia[1]-1;iy<=ia[1]+1;iy++){
        if(iy<0 || iy>=numOfVoxel[1]) continue;
        for(int iz=ia[2]-1;iz<=ia[2]+1;iz++){
          if(iz<0 || iz>=numOfVoxel[2]) continue;

          for(int it=0;it<numOfNodeInBox[ix][iy][iz][ibeam2];it++){

            ic2 = nodeInBox[ix][iy][iz][ibeam2][it];

            ContactInBox(beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component,ic1,ic2,beam[ibeam1],beam[ibeam2],beamContactHistory[ibeam1][ibeam2-ibeam1][ic1].contactHistoryMAP,kcon,E_astalisk);
          }
        }
      }
    }
  }

  //assembly
  int icc;
  for(int ic1=0;ic1<beam[ibeam1].numOfElm;ic1++){
    for(int ic2=0;ic2<beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component.size();ic2++){
      icc=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].ic;
      for(int j=0;j<3;j++){
        beam[ibeam1].fcon[beam[ibeam1].ie[ic1][0]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].fn[0][j];
        beam[ibeam1].fcon[beam[ibeam1].ie[ic1][1]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].fn[1][j];
        beam[ibeam2].fcon[beam[ibeam2].ie[icc][0]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].fn[2][j];
        beam[ibeam2].fcon[beam[ibeam2].ie[icc][1]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].fn[3][j];
        beam[ibeam1].fcon[beam[ibeam1].ie[ic1][0]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].ft[0][j];
        beam[ibeam1].fcon[beam[ibeam1].ie[ic1][1]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].ft[1][j];
        beam[ibeam2].fcon[beam[ibeam2].ie[icc][0]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].ft[2][j];
        beam[ibeam2].fcon[beam[ibeam2].ie[icc][1]][j]+=beamContactTMP[ibeam1][ibeam2-ibeam1][ic1].component[ic2].ft[3][j];
      }
    }
  }

}

// #########################################################
/**
 * @brief contact force
 * @param [out]    component    Temporal data strage for contact
 * @param [in]     ic1          element number
 * @param [in]     ic2          element number
 * @param [in]     beam1        target beam1 (beam[ibeam1])
 * @param [in]     beam2        target beam2 (beam[ibeam2])
 * @param [in]     map          contact history <contact element, frictional displacement[3]>
 * @param [in]     kcon         elastic contact coefficient
 * @param [in]     E_astalisk   pseudo-elastic coefficient (please see Heltz contact theory)
 */
void multipleBeamSimulator::ContactInBox(std::vector<BeamContactComponent> &component,const int ic1,const int ic2,coRotationalBeam &beam1,coRotationalBeam &beam2,
      std::unordered_map<int,ContactDisplacement> map,const double kcon,const double E_astalisk)
{
  double distance,s,t,normal[3];
  double collision=1e-2;
  double va1[3],va2[3],vn[3],us[3],fn[4][3],ft[4][3];
  BeamContactComponent BCTmp;

  distance=beam2beam_contact_search(beam1.x_p,beam2.x_p,beam1,beam2,ic1,ic2,s,t,normal);

  if(distance>beam1.rad+beam2.rad*1.1e0) return;
  if(distance>beam1.rad+beam2.rad && map.count(ic2)==0) return;

  double M1=(1e0-s)*beam1.Ma[beam1.ie[ic1][0]]+s*beam1.Ma[beam1.ie[ic1][1]];
  double M2=(1e0-t)*beam2.Ma[beam2.ie[ic2][0]]+t*beam2.Ma[beam2.ie[ic2][1]];

  double vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(5e-1*(M1+M2)*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));

  BCTmp.ic=ic2;

  if(map.count(ic2)!=0){
    for(int j=0;j<3;j++) us[j]=map[ic2].u_static[j];
    BCTmp.stickSlip=true;
  }else{
    BCTmp.stickSlip=false;
  }

  for(int j=0;j<3;j++) normal[j]=normal[j]/distance;
  for(int j=0;j<3;j++){
    va1[j]=(1e0-s)*beam1.v_p[beam1.ie[ic1][0]][j]+s*beam1.v_p[beam1.ie[ic1][1]][j];
    va2[j]=(1e0-t)*beam2.v_p[beam2.ie[ic2][0]][j]+t*beam2.v_p[beam2.ie[ic2][1]][j];
    vn[j]=va2[j]-va1[j];
  }

  if(distance<beam1.rad+beam2.rad){
    calc_beamContactForce(fn,beam1,beam2,s,t,kcon,vcon,distance,normal,vn);
    calc_beamContactForce_friction(ft,beam1,beam2,s,t,kcon,distance,normal,vn,us,E_astalisk,BCTmp.stickSlip);
    for(int i=0;i<4;i++){
      for(int j=0;j<3;j++){
        BCTmp.fn[i][j]=fn[i][j];
        BCTmp.ft[i][j]=ft[i][j];
      }
    }
    for(int j=0;j<3;j++) BCTmp.u_static[j]=us[j];
    component.push_back(BCTmp);

  }else if(map.count(ic2)!=0){
    calc_beamContactForce_friction(ft,beam1,beam2,s,t,kcon,distance,normal,vn,us,E_astalisk,BCTmp.stickSlip);
    for(int i=0;i<4;i++){
      for(int j=0;j<3;j++){
        BCTmp.fn[i][j]=0e0;
        BCTmp.ft[i][j]=ft[i][j];
      }
    }
    for(int j=0;j<3;j++) BCTmp.u_static[j]=us[j];
    component.push_back(BCTmp);
  }
}

// #########################################################
/**
 * @brief contact force (normal direction)
 * @param [out]    fn           Normal contact force
 * @param [in]     beam1        target beam1 (beam[ibeam1])
 * @param [in]     beam2        target beam2 (beam[ibeam2])
 * @param [in]     s            shape function (1st-order) of beam1
 * @param [in]     t            shape function (1st-order) of beam2
 * @param [in]     kcon         elastic contact coefficient
 * @param [in]     vcon         viscos contact coefficient
 * @param [in]     distance     distance between beam1 and bram2
 * @param [in]     normal       unit normal vector
 * @param [in]     vn           relative velocity vector along normal direction
 */
void multipleBeamSimulator::calc_beamContactForce(double (&fn)[4][3],const coRotationalBeam &beam1,const coRotationalBeam &beam2,const double s,const double t,
  const double kcon,const double vcon,const double distance,const double (&normal)[3],const double (&vn)[3])
{
  double deflection=fabs(beam1.rad+beam2.rad-distance);
  double fnormal=kcon*deflection*5e0;

  double vtmp=vn[0]*normal[0]+vn[1]*normal[1]+vn[2]*normal[2];
  double vnormal=vcon*vtmp;

  for(int j=0;j<3;j++){
    fn[0][j]=-(fnormal-vnormal)*normal[j]*(1e0-s);
    fn[1][j]=-(fnormal-vnormal)*normal[j]*s;
    fn[2][j]=(fnormal-vnormal)*normal[j]*(1e0-t);
    fn[3][j]=(fnormal-vnormal)*normal[j]*t;
  }

}

// #########################################################
/**
 * @brief contact force (tangential or frictional)
 * @param [out]    ft           Tangential (frictional) contact force
 * @param [in]     beam1        target beam1 (beam[ibeam1])
 * @param [in]     beam2        target beam2 (beam[ibeam2])
 * @param [in]     s            shape function (1st-order) of beam1
 * @param [in]     t            shape function (1st-order) of beam2
 * @param [in]     kcon         elastic contact coefficient
 * @param [in]     distance     distance between beam1 and bram2
 * @param [in]     normal       unit normal vector
 * @param [in]     vn           relative velocity vector along normal direction
 * @param [out]    us           frictional displacement (only acitive in stick mode)
 * @param [in]     E_astalisk   pseudo-elastic coefficient (please see Heltz contact theory)
 * @param [out]    stickSlip    stick/slip mode info.
 */
void multipleBeamSimulator::calc_beamContactForce_friction(double (&ft)[4][3],const coRotationalBeam &beam1,const coRotationalBeam &beam2,const double s,const double t,const double kcon,const double distance,
const double (&normal)[3],const double (&vn)[3],double (&us)[3],const double E_astalisk,bool &stickSlip)
{
  double mu_slip=9e-1*mu_static_coil; //Vetter, 2015
  double k_static=E_astalisk*(beam1.rad+beam2.rad)*5e-1;
  double veps=1e-8*minimumSlidingVelocity;
  double vt[3],ftmp,ft_tmp[3];

  double vtmp=vn[0]*normal[0]+vn[1]*normal[1]+vn[2]*normal[2];

  for(int j=0;j<3;j++) vt[j]=vn[j]-vtmp*normal[j];
  double vt_abs=sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);

  double deflection=fabs(beam1.rad+beam2.rad-distance);
  double fnormal=kcon*deflection;
  double ft_abs=mu_slip*fnormal;

  if(vt_abs>=minimumSlidingVelocity){
    for(int j=0;j<3;j++) ft_tmp[j]=-ft_abs*vt[j]/(vt_abs+veps);
    stickSlip=false;
  }else{
    if(stickSlip==false){   //initialize
      for(int j=0;j<3;j++) us[j]=(ft_abs/(vt_abs+veps)/k_static-stickFrictionViscosCoefficient)*vt[j];
    }else{
      for(int j=0;j<3;j++) us[j]+=vt[j]*dt;
    }
    for(int j=0;j<3;j++) ft_tmp[j]=-k_static*us[j]-stickFrictionViscosCoefficient*k_static*vt[j];
    ftmp=sqrt(ft_tmp[0]*ft_tmp[0]+ft_tmp[1]*ft_tmp[1]+ft_tmp[2]*ft_tmp[2]);
    stickSlip=true;
    if(ftmp>mu_static_coil*fnormal){
      for(int j=0;j<3;j++) ft_tmp[j]=-ft_abs*vt[j]/(vt_abs+veps);
      stickSlip=false;
    }
  }
  for(int j=0;j<3;j++){
    ft[0][j]=-ft_tmp[j]*(1e0-s);
    ft[1][j]=-ft_tmp[j]*s;
    ft[2][j]=ft_tmp[j]*(1e0-t);
    ft[3][j]=ft_tmp[j]*t;
  }
}

// #########################################################
/**
 * @brief contact search (Wriggers & Zabrize? 2000)
 * @param [in]     x1           Tangential (frictional) contact force
 * @param [in]     x2           target beam1 (beam[ibeam1])
 * @param [in]     beam1        target beam1 (beam[ibeam1])
 * @param [in]     beam2        target beam2 (beam[ibeam2])
 * @param [in]     ic1          element number
 * @param [in]     ic2          element number
 * @param [out]    s            shape function (1st-order) of beam1
 * @param [out]    t            shape function (1st-order) of beam2
 * @param [out]    g            normal vector
 */
double multipleBeamSimulator::beam2beam_contact_search(const DOUBLEARRAY2 &x1,const DOUBLEARRAY2 &x2,
      coRotationalBeam &beam1,coRotationalBeam &beam2,const int ic1,const int ic2,double &s,double &t,double (&g)[3])
 {
  double eps=1e-15;
  double a[3],b[3],c[3];

  for(int j=0;j<3;j++){
    a[j]=x1[beam1.ie[ic1][1]][j]-x1[beam1.ie[ic1][0]][j];
    b[j]=x2[beam2.ie[ic2][1]][j]-x2[beam2.ie[ic2][0]][j];
    c[j]=x1[beam1.ie[ic1][0]][j]-x2[beam2.ie[ic2][0]][j];
  }

  double aa=0e0,bb=0e0,ab=0e0,ac=0e0,bc=0e0;
  for(int j=0;j<3;j++){
    aa+=a[j]*a[j];
    ab+=a[j]*b[j];
    bb+=b[j]*b[j];
    ac+=a[j]*c[j];
    bc+=b[j]*c[j];
  }

  double tmp=ab*ab-aa*bb;

  if(fabs(tmp)<eps){ //if parallel
    s=5e-1;
    t=5e-1;
  }else{
    s=1e0/tmp*(bb*ac-ab*bc);
    t=1e0/tmp*(ab*ac-aa*bc);
    if(s<0e0) s=0e0;
    if(1e0<s) s=1e0;
    if(t<0e0) t=0e0;
    if(1e0<t) t=1e0;
  }

  double xa1[3],xa2[3];
  for(int j=0;j<3;j++){
    xa1[j]=(1e0-s)*x1[beam1.ie[ic1][0]][j]+s*x1[beam1.ie[ic1][1]][j];
    xa2[j]=(1e0-t)*x2[beam2.ie[ic2][0]][j]+t*x2[beam2.ie[ic2][1]][j];
    g[j]=xa2[j]-xa1[j];
  }
  return sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
 }



// #########################################################
/**
 * @brief contact search contact beam
 */
 // void multipleBeamSimulator::search_contact_beam(const int ibeam1, const int ibeam2)
 // {
 //  double tmp0=1e2;
 //  int tmp1, tmp2;
 //  double tmp3, tmp4;
 //  double distance,s,t,normal[3];

 //  for(int ic1=0; ic1<beam[ibeam1].numOfElm; ic1++){
 //    for(int j=0; j<2; j++){
 //      contact_beam[ibeam1][ibeam2][ic1][j] = -1;
 //    }
 //  }


 //  for(int ic1=0; ic1<beam[ibeam1].numOfElm; ic1++){
 //    tmp0=1e2;
 //    for(int ic2=0; ic2<beam[ibeam2].numOfElm; ic2++){
 //      distance = beam2beam_contact_search(beam[ibeam1].x, beam[ibeam2].x, beam[ibeam1],beam[ibeam2],ic1,ic2,s,t,normal);
 //      if(distance < tmp0){
 //        tmp0 = distance;
 //        tmp1 = ic1;
 //        tmp2 = ic2;
 //        tmp3 = s;
 //        tmp4 = t;
 //      }
 //    }

 //    if(tmp0<1e-2){
 //      if(tmp3<1e-10 && tmp4<1e-10){
 //        contact_beam[ibeam1][ibeam2][ic1][0] = tmp1;
 //        contact_beam[ibeam1][ibeam2][ic1][1] = tmp2;
 //        // printf("ibeam1=%2d ibeam2=%2d ic1=%3d ic2=%3d s=%e t=%e distance=%e  contact_beam[0]=%d  contact_beam[1]=%d\n", 
 //        //          ibeam1,ibeam2,tmp1,tmp2,tmp3,tmp4,tmp0,contact_beam[ibeam1][ibeam2][ic1][0],contact_beam[ibeam1][ibeam2][ic1][1]);
 //      }

 //      if(tmp3>1e-10 && tmp3<1e0){
 //        if(tmp4>1e-10 && tmp4<1e0){
 //          contact_beam[ibeam1][ibeam2][ic1][0] = tmp1;
 //          contact_beam[ibeam1][ibeam2][ic1][1] = tmp2;
 //          // printf("ibeam1=%2d ibeam2=%2d ic1=%3d ic2=%3d s=%e t=%e distance=%e  contact_beam[0]=%d  contact_beam[1]=%d\n", 
 //          //          ibeam1,ibeam2,tmp1,tmp2,tmp3,tmp4,tmp0,contact_beam[ibeam1][ibeam2][ic1][0],contact_beam[ibeam1][ibeam2][ic1][1]);
 //        }
 //      }      
 //    }
 //  }
 // }


 // #########################################################
/**
 * @brief contact search contact beam
 */
//  void multipleBeamSimulator::search_contact_beam_2(const int ibeam1, const int ic1)
//  {
//   double tmp0=1e2;
//   int tmp1, tmp2;
//   double tmp3, tmp4;
//   int tmp5;
//   double distance,s,t,normal[3];

//   for(int j=0; j<3; j++){
//     beam[ibeam1].contact_beam2[ic1][j] = -1;
//   }

//   tmp0 = 1e2;
//   for(int ibeam2=NW1; ibeam2<NW1+NW2; ibeam2++){
//     for(int ic2=0; ic2<beam[ibeam2].numOfElm; ic2++){
//       distance = beam2beam_contact_search(beam[ibeam1].x, beam[ibeam2].x, beam[ibeam1],beam[ibeam2],ic1,ic2,s,t,normal);
//       if(distance < tmp0){
//         tmp0 = distance;
//         tmp1 = ic1;
//         tmp2 = ic2;
//         tmp3 = s;
//         tmp4 = t;
//         tmp5 = ibeam2;
//       }
//     }
//   }

//   if(tmp0<1e-2){
//     if(tmp3<1e-10 && tmp4<1e-10){
//       beam[ibeam1].contact_beam2[ic1][0] = tmp1;
//       beam[ibeam1].contact_beam2[ic1][1] = tmp2;
//       beam[ibeam1].contact_beam2[ic1][2] = tmp5;
//        //printf("ibeam1=%2d ibeam2=%2d ic1=%3d ic2=%3d s=%8e t=%8e distance=%e  contact_beam[0]=%d  contact_beam[1]=%d  contact_beam[2]=%d\n", 
//         //        ibeam1,tmp5,tmp1,tmp2,tmp3,tmp4,tmp0,beam[ibeam1].contact_beam2[ic1][0],beam[ibeam1].contact_beam2[ic1][1], beam[ibeam1].contact_beam2[ic1][2]);
//     }

//     if(tmp3>1e-10 && tmp3<1e0){
//       if(tmp4>1e-10 && tmp4<1e0){
//         beam[ibeam1].contact_beam2[ic1][0] = tmp1;
//         beam[ibeam1].contact_beam2[ic1][1] = tmp2;
//         beam[ibeam1].contact_beam2[ic1][2] = tmp5;
//          //printf("ibeam1=%2d ibeam2=%2d ic1=%3d ic2=%3d s=%8e t=%8e distance=%e  contact_beam[0]=%d  contact_beam[1]=%d  contact_beam[2]=%d\n", 
//            //       ibeam1,tmp5,tmp1,tmp2,tmp3,tmp4,tmp0,beam[ibeam1].contact_beam2[ic1][0],beam[ibeam1].contact_beam2[ic1][1], beam[ibeam1].contact_beam2[ic1][2]);
//       }
//     }      
//   }
//  }

  // #########################################################
/**
 * @brief contact search contact beam
 */
//  void multipleBeamSimulator::wire_cross_force(const int ibeam1, const int ibeam2, const int ic1, const int ic2, const double exf_norm, const DOUBLEARRAY1 wire_cross)
//  {
//   double distance,s,t,normal[3];
//   double er[3];
//   double er_norm;
//   double u_norm = 1e0;

//   distance = beam2beam_contact_search(beam[ibeam1].x, beam[ibeam2].x, beam[ibeam1],beam[ibeam2],ic1,ic2,s,t,normal);

//   er[0] = beam[ibeam1].x[ic1][0];
//   er[1] = beam[ibeam1].x[ic1][1];
//   er[2] = 0e0;
//   er_norm = sqrt(er[0]*er[0] + er[1]*er[1]);
//   for(int j=0; j<3; j++){
//     beam[ibeam1].v[beam[ibeam1].ie[ic1][0]][j] = (1e0-s) * u_norm * wire_cross[ic1] * er[j];
//     beam[ibeam1].v[beam[ibeam1].ie[ic1][1]][j] = s * u_norm * wire_cross[ic1] * er[j];

//     beam[ibeam2].v[beam[ibeam2].ie[ic2][0]][j] = (1e0-t) * u_norm * wire_cross[ic1]* (-1e0) * er[j];
//     beam[ibeam2].v[beam[ibeam2].ie[ic2][1]][j] = t * u_norm * wire_cross[ic1] * (-1e0) * er[j];
//   }  
//  }


// #########################################################
/**
 * @brief contact force
 */
// void multipleBeamSimulator::selfContact_APC(const int &ibeam)
// {
//   double distance,s,t;

//   for(int ic1=0;ic1<beam[ibeam].numOfElm;ic1++){
//     for(int ic2=ic1+2;ic2<beam[ibeam].numOfElm;ic2++){
//       distance=beam2beam_contact_search(beam[ibeam].x_p,beam[ibeam].x_p,beam[ibeam],beam[ibeam],ic1,ic2,s,t);
//       if(distance<beam[ibeam].rad+beam[ibeam].rad){
//         calc_beamSelfContactForce(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,ic1,ic2,s,t);
//       }
//     }
//   }
// }

// #########################################################
/**
 * @brief contact force
 */
// void multipleBeamSimulator::selfContact(const int &ibeam)
// {
//   double distance,s,t,norma;[3];

//   for(int ic1=0;ic1<beam[ibeam].numOfElm;ic1++){
//     for(int ic2=ic1+2;ic2<beam[ibeam].numOfElm;ic2++){
//       distance=beam2beam_contact_search(beam[ibeam].x,beam[ibeam].x,beam[ibeam],beam[ibeam],ic1,ic2,s,t);
//       if(distance<beam[ibeam].rad+beam[ibeam].rad){
//         calc_beamSelfContactForce(beam[ibeam],beam[ibeam].x,beam[ibeam].v,ic1,ic2,s,t);
//       }
//     }
//   }
// }
// #########################################################
/**
 * @brief contact force
 */
// void multipleBeamSimulator::contact_beams_APC(const int &ibeam1,const int &ibeam2)
// {
//   int check;
//   double distance,s,t;

//   for(int ic1=0;ic1<beam[ibeam1].numOfElm;ic1++){

//     check=0;

//     for(int ic2=0;ic2<beam[ibeam2].numOfElm;ic2++){
//       distance=beam2beam_contact_search(beam[ibeam1].x_p,beam[ibeam2].x_p,beam[ibeam1],beam[ibeam2],ic1,ic2,s,t);
//       if(distance<beam[ibeam1].rad+beam[ibeam2].rad){

//         //to avoid duplicated contact detection
//         if(t<1e-15 && check==1){
//           //printf("beam=%d beam=%d elm=%d elm-%d\n",ibeam1,ibeam2,ic1,ic2);
//           continue;
//         }
//         if(abs(t-1e0)<1e-15){
//           check=1;
//         }else{
//           check=0;
//         }

//         calc_beam2beamContactForce(beam[ibeam1],beam[ibeam2],beam[ibeam1].x_p,beam[ibeam2].x_p,beam[ibeam1].v_p,beam[ibeam2].v_p,ic1,ic2,s,t);
//       }
//     }
//   }
// }

// #########################################################
/**
 * @brief contact force
 */
// void multipleBeamSimulator::contact_beams(const int &ibeam1,const int &ibeam2)
// {
//   double distance,s,t;

//   for(int ic1=0;ic1<beam[ibeam1].numOfElm;ic1++){
//     for(int ic2=0;ic2<beam[ibeam2].numOfElm;ic2++){
//       distance=beam2beam_contact_search(beam[ibeam1].x,beam[ibeam2].x,beam[ibeam1],beam[ibeam2],ic1,ic2,s,t);
//       if(distance<beam[ibeam1].rad+beam[ibeam2].rad){
//         calc_beam2beamContactForce(beam[ibeam1],beam[ibeam2],beam[ibeam1].x,beam[ibeam2].x,beam[ibeam1].v,beam[ibeam2].v,ic1,ic2,s,t);
//       }
//     }
//   }
// }
