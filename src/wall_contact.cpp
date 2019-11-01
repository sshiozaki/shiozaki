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
 * @file wall_contact.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

// #########################################################
/**
 * @brief contact force
 * @param [in]     ibeam    target beam numberd
 */
void multipleBeamSimulator::closedSurfaceContact_APC(const int ibeam)
{
  double PI=acos(-1.0e0);
  int ix[3];
  double g[3],length,normal[3];
  double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
        +(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0)*1e1);

 //double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0));

  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;
  double vcon;

  #pragma omp parallel for private(length,g,normal,ix,vcon)
  for(int ic=0;ic<beam[ibeam].numOfNode;ic++){
    if(beam[ibeam].mask1[ic][1]==0) continue;
    length=contactDetection_SDF(g,ix,beam[ibeam],beam[ibeam].x_p,ic);
    if(length-beam[ibeam].rad>beam[ibeam].rad*1e-1){
      wallContact[ibeam][ic]=0;
      continue;
    }
    vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(beam[ibeam].Ma[ic]*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));
    calc_unitNormalVector(normal,g,ix);

    if(length<beam[ibeam].rad){
      calc_wallContactForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,kcon,vcon,normal,length,ic);
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }else if(wallContact[ibeam][ic]==1){
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }

    for(int j=0; j<3; j++){
      beam[ibeam].f_wall[ic][j] = beam[ibeam].fcon[ic][j];
      beam[ibeam].ft_wall[ic][j] = beam[ibeam].f_wall[ic][j] - beam[ibeam].fn_wall[ic][j];
    }
  }

}

// #########################################################
/**
 * @brief contact force
 * @param [in]     ibeam    target beam numberd
 */
void multipleBeamSimulator::closedSurfaceContact_APC_catheter(const int ibeam)
{
  double PI=acos(-1.0e0);
  int ix[3];
  double g[3],length,normal[3];

  //double EA_catheter=30e3*PI*pow(beam[ibeam].rad,2e0);
  double EA_catheter=30e3;
  double Poisson_catheter=0.3;

  double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
        +(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0)*1e1);

  // double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
  //       +(1e0-Poisson_catheter*Poisson_catheter)/(EA_catheter*PI*pow(beam[ibeam].rad,2e0));

  //printf("kcon= %e kcon2= %e Poisson= %e %e %e\n", kcon, kcon2, beam[0].Poisson, beam[0].EA/(PI*pow(beam[0].rad,2e0)), beam[0].EA);

  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;

  double vcon;

  #pragma omp parallel for private(length,g,normal,ix,vcon)
  for(int ic=0;ic<beam[ibeam].numOfNode;ic++){
    if(beam[ibeam].mask1[ic][1]==0) continue;
    length=contactDetection_SDF_catheter(g,ix,beam[ibeam],beam[ibeam].x_p,ic);
    if(length-beam[ibeam].rad>beam[ibeam].rad*1e-1){
      wallContact[ibeam][ic]=0;
      continue;
    }
    vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(beam[ibeam].Ma[ic]*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));
    calc_unitNormalVector_catheter(normal,g,ix);

    if(length<beam[ibeam].rad){
      calc_wallContactForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,kcon,vcon,normal,length,ic);
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }else if(wallContact[ibeam][ic]==1){
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }
  }

}

// #########################################################
/**
 * @brief contact force
 * @param [in]     ibeam    target beam numberd
 */
void multipleBeamSimulator::closedSurfaceContact_APC_2(const int ibeam, const int release_node_number)
{
  double PI=acos(-1.0e0);
  int ix[3];
  double g[3],length,normal[3];
  double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
        +(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0)*1e1);

  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;
  double vcon;

  #pragma omp parallel for private(length,g,normal,ix,vcon)
  for(int ic=beam[ibeam].numOfNode-1; ic>beam[ibeam].numOfNode-1-release_node_number ;ic--){           //caution!!!!!!!!!!!
    if(beam[ibeam].mask1[ic][1]==0) continue;
    length=contactDetection_SDF(g,ix,beam[ibeam],beam[ibeam].x_p,ic);
    if(length-beam[ibeam].rad>beam[ibeam].rad*1e-1){
      wallContact[ibeam][ic]=0;
      continue;
    }
    vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(beam[ibeam].Ma[ic]*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));
    calc_unitNormalVector(normal,g,ix);

    if(length<beam[ibeam].rad){
      calc_wallContactForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,kcon,vcon,normal,length,ic);
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }else if(wallContact[ibeam][ic]==1){
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }

    for(int j=0; j<3; j++){
      beam[ibeam].f_wall[ic][j] = beam[ibeam].fcon[ic][j];
      beam[ibeam].ft_wall[ic][j] = beam[ibeam].f_wall[ic][j] - beam[ibeam].fn_wall[ic][j];
    }
  }

}

// #########################################################
/**
 * @brief contact force
 * @param [in]     ibeam    target beam numberd
 */
void multipleBeamSimulator::closedSurfaceContact_APC_catheter_2(const int ibeam, const int release_node_number)
{
  double PI=acos(-1.0e0);
  int ix[3];
  double g[3],length,normal[3];

  double EA_catheter=30e3;
  double Poisson_catheter=0.3;

  double kcon=(1e0-beam[ibeam].Poisson*beam[ibeam].Poisson)/(beam[ibeam].EA*PI*pow(beam[ibeam].rad,2e0))
        +(1e0-Poisson_catheter*Poisson_catheter)/(EA_catheter*PI*pow(beam[ibeam].rad,2e0));

  double E_astalisk=1e0/kcon;
  kcon=2.5e-1*acos(-1.0e0)*E_astalisk;

  double vcon;

  #pragma omp parallel for private(length,g,normal,ix,vcon)
  for(int ic=0;ic<beam[ibeam].numOfNode-release_node_number;ic++){
    if(beam[ibeam].mask1[ic][1]==0) continue;
    length=contactDetection_SDF_catheter(g,ix,beam[ibeam],beam[ibeam].x_p,ic);
    if(length-beam[ibeam].rad>beam[ibeam].rad*1e-1){
      wallContact[ibeam][ic]=0;
      continue;
    }
    vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(beam[ibeam].Ma[ic]*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));
    calc_unitNormalVector_catheter(normal,g,ix);

    if(length<beam[ibeam].rad){
      calc_wallContactForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,kcon,vcon,normal,length,ic);
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }else if(wallContact[ibeam][ic]==1){
      calc_wallFrictionForce_SDF(beam[ibeam],beam[ibeam].x_p,beam[ibeam].v_p,E_astalisk,kcon,normal,length,ibeam,ic);
    }

    for(int j=0; j<3; j++){
      beam[ibeam].f_wall[ic][j] = beam[ibeam].fcon[ic][j];
      beam[ibeam].ft_wall[ic][j] = beam[ibeam].f_wall[ic][j] - beam[ibeam].fn_wall[ic][j];
    }
  }

}

// #########################################################
/**
 * @brief calc unit normal vector from SDF function
  * @param [out]    normal    target beam number
  * @param [in]     g         position vector in voxel
  * @param [in]     ix        voxel number [x,y,z]
 */
void multipleBeamSimulator::calc_unitNormalVector(double (&normal)[3],const double (&g)[3],const int (&ix)[3])
{
  double dNdr[8][3],tmp;

  wall[0].C3D8_dNdr(dNdr,g[0],g[1],g[2]);
  for(int i=0;i<3;++i){
    normal[i] = dNdr[0][i] * wall[0].SDF[ix[0]][ix[1]][ix[2]]    + dNdr[1][i] * wall[0].SDF[ix[0]+1][ix[1]][ix[2]]
          + dNdr[2][i] * wall[0].SDF[ix[0]+1][ix[1]+1][ix[2]]  + dNdr[3][i] * wall[0].SDF[ix[0]][ix[1]+1][ix[2]]
          + dNdr[4][i] * wall[0].SDF[ix[0]][ix[1]][ix[2]+1]   + dNdr[5][i] * wall[0].SDF[ix[0]+1][ix[1]][ix[2]+1]
          + dNdr[6][i] * wall[0].SDF[ix[0]+1][ix[1]+1][ix[2]+1] + dNdr[7][i] * wall[0].SDF[ix[0]][ix[1]+1][ix[2]+1];
  }
  tmp = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  for(int i=0;i<3;i++) normal[i] = -normal[i]/tmp;
}

// #########################################################
/**
 * @brief calc unit normal vector from SDF function
  * @param [out]    normal    target beam number
  * @param [in]     g         position vector in voxel
  * @param [in]     ix        voxel number [x,y,z]
 */
void multipleBeamSimulator::calc_unitNormalVector_catheter(double (&normal)[3],const double (&g)[3],const int (&ix)[3])
{
  double dNdr[8][3],tmp;

  wall[1].C3D8_dNdr(dNdr,g[0],g[1],g[2]);
  for(int i=0;i<3;++i){
    normal[i] = dNdr[0][i] * wall[1].SDF[ix[0]][ix[1]][ix[2]]    + dNdr[1][i] * wall[1].SDF[ix[0]+1][ix[1]][ix[2]]
          + dNdr[2][i] * wall[1].SDF[ix[0]+1][ix[1]+1][ix[2]]  + dNdr[3][i] * wall[1].SDF[ix[0]][ix[1]+1][ix[2]]
          + dNdr[4][i] * wall[1].SDF[ix[0]][ix[1]][ix[2]+1]   + dNdr[5][i] * wall[1].SDF[ix[0]+1][ix[1]][ix[2]+1]
          + dNdr[6][i] * wall[1].SDF[ix[0]+1][ix[1]+1][ix[2]+1] + dNdr[7][i] * wall[1].SDF[ix[0]][ix[1]+1][ix[2]+1];
  }
  tmp = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  for(int i=0;i<3;i++) normal[i] = -normal[i]/tmp;
}

// #########################################################
/**
 * @brief contact detection from SDF function
  * @param [in]     g         position vector in voxel
  * @param [in]     ix        voxel number [x,y,z]
  * @param [in]     beam      target beam
  * @param [in]     x         node position vector
  * @param [in]     ic        node number
 */
double multipleBeamSimulator::contactDetection_SDF(double (&g)[3],int (&ix)[3],coRotationalBeam &beam,const DOUBLEARRAY2 &x,const int ic)
{
  double N[8];

  for(int j=0;j<3;j++){
    g[j] = (beam.x[ic][j]-wall[0].origin[j])/wall[0].dx[j];
    ix[j] = (int)g[j];
    g[j] = 2e0*(g[j]-(double)ix[j])-1e0;
  }

  wall[0].C3D8_N(N,g[0],g[1],g[2]);

  double length = N[0] * wall[0].SDF[ix[0]][ix[1]][ix[2]]    + N[1] * wall[0].SDF[ix[0]+1][ix[1]][ix[2]]
          + N[2] * wall[0].SDF[ix[0]+1][ix[1]+1][ix[2]]  + N[3] * wall[0].SDF[ix[0]][ix[1]+1][ix[2]]
          + N[4] * wall[0].SDF[ix[0]][ix[1]][ix[2]+1]   + N[5] * wall[0].SDF[ix[0]+1][ix[1]][ix[2]+1]
          + N[6] * wall[0].SDF[ix[0]+1][ix[1]+1][ix[2]+1] + N[7] * wall[0].SDF[ix[0]][ix[1]+1][ix[2]+1];
  return length;
}

// #########################################################
/**
 * @brief contact detection from SDF function
  * @param [in]     g         position vector in voxel
  * @param [in]     ix        voxel number [x,y,z]
  * @param [in]     beam      target beam
  * @param [in]     x         node position vector
  * @param [in]     ic        node number
 */
double multipleBeamSimulator::contactDetection_SDF_catheter(double (&g)[3],int (&ix)[3],coRotationalBeam &beam,const DOUBLEARRAY2 &x,const int ic)
{
  double N[8];

  for(int j=0;j<3;j++){
    g[j] = (beam.x[ic][j]-wall[1].origin[j])/wall[1].dx[j];
    ix[j] = (int)g[j];
    g[j] = 2e0*(g[j]-(double)ix[j])-1e0;
  }

  wall[1].C3D8_N(N,g[0],g[1],g[2]);

  double length = N[0] * wall[1].SDF[ix[0]][ix[1]][ix[2]]    + N[1] * wall[1].SDF[ix[0]+1][ix[1]][ix[2]]
          + N[2] * wall[1].SDF[ix[0]+1][ix[1]+1][ix[2]]  + N[3] * wall[1].SDF[ix[0]][ix[1]+1][ix[2]]
          + N[4] * wall[1].SDF[ix[0]][ix[1]][ix[2]+1]   + N[5] * wall[1].SDF[ix[0]+1][ix[1]][ix[2]+1]
          + N[6] * wall[1].SDF[ix[0]+1][ix[1]+1][ix[2]+1] + N[7] * wall[1].SDF[ix[0]][ix[1]+1][ix[2]+1];
  return length;
}
// #########################################################
/**
 * @brief contact normal force from SDF function
  * @param [inout]  beam      target beam
  * @param [in]     x         node position vector
  * @param [in]     v         node velocity vector
  * @param [in]     kcon      elastic contact coefficient
  * @param [in]     vcon      viscos contact coefficient
  * @param [in]     normal    unit normal vector
  * @param [in]     length    distance between wall and beam node
  * @param [in]     ic        node number
 */
void multipleBeamSimulator::calc_wallContactForce_SDF(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,const double kcon,const double vcon,
                            const double (&normal)[3],const double length,const int ic)
{
  //normal force
  double delta = fabs(beam.rad-length);
  double fnormal=kcon*delta;

  double vtmp = beam.v[ic][0]*normal[0] + beam.v[ic][1]*normal[1] + beam.v[ic][2]*normal[2];
  double vnormal=vcon*vtmp;

  for(int j=0;j<3;j++){
    beam.fcon[ic][j]-=fnormal*normal[j];
    beam.fcon[ic][j]-=vnormal*normal[j];
  }

  for(int j=0; j<3; j++){
    beam.fn_wall[ic][j] = beam.fcon[ic][j];
  }
}
// #########################################################
/**
 * @brief contact force
  * @param [inout]  beam        target beam
  * @param [in]     x           node position vector
  * @param [in]     v           node velocity vector
 * @param [in]     E_astalisk   pseudo-elastic coefficient (please see Heltz contact theory)
  * @param [in]     kcon        elastic contact coefficient
  * @param [in]     normal      unit normal vector
  * @param [in]     length      distance between wall and beam node
  * @param [in]     ibeam       target beam number
  * @param [in]     ic          node number
 */
void multipleBeamSimulator::calc_wallFrictionForce_SDF(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,
                    const double E_astalisk,const double kcon,const double (&normal)[3],const double length,const int ibeam,const int ic)
 {
  double mu_slip=9e-1*mu_static_aneurysm;
  double k_static=E_astalisk*beam.rad;
  double veps=1e-8*minimumSlidingVelocity;
  double vn[3],ut,vt[3],ft[3],ftmp;

  double vtmp = beam.v[ic][0]*normal[0] + beam.v[ic][1]*normal[1] + beam.v[ic][2]*normal[2];

  //normal force
  double delta = fabs(beam.rad-length);
  double fnormal=kcon*delta;

  for(int j=0;j<3;j++) vt[j]=beam.v[ic][j]-vtmp*normal[j];
  double vt_abs=sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
  double ft_abs=mu_slip*fnormal;

  if(vt_abs>=minimumSlidingVelocity){
    for(int j=0;j<3;j++) ft[j]=-ft_abs*vt[j]/(vt_abs+veps);
    wallContact[ibeam][ic]=0;
  }else{
    if(wallContact[ibeam][ic]==0){  //initialize
      for(int j=0;j<3;j++) u_stat_wall[ibeam][ic][j]=(ft_abs/(vt_abs+veps)/k_static-stickFrictionViscosCoefficient)*vt[j];
    }else{
      for(int j=0;j<3;j++) u_stat_wall_tmp[ibeam][ic][j]=u_stat_wall[ibeam][ic][j]+vt[j]*dt;
    }

    for(int j=0;j<3;j++) ft[j]=-k_static*u_stat_wall_tmp[ibeam][ic][j]-stickFrictionViscosCoefficient*k_static*vt[j];
    ftmp=sqrt(ft[0]*ft[0]+ft[1]*ft[1]+ft[2]*ft[2]);
    if(ftmp>mu_static_aneurysm*fnormal){
      for(int j=0;j<3;j++) ft[j]=-ft_abs*vt[j]/(vt_abs+veps);
      wallContact[ibeam][ic]=0;
    }else{
      wallContact[ibeam][ic]=1;
    }
  }

  for(int j=0;j<3;j++) beam.fcon[ic][j]+=ft[j];

 }

// #########################################################
/**
 * @brief wall contact force
  * @param [inout]  beam        target beam
  * @param [in]     x           node position vector
  * @param [in]     v           node velocity vector
  * @param [in]     ic          node number
  * @param [in]     Radius      beam radius
 */
void multipleBeamSimulator::calc_wallContactForce(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,const int ic,const double Radius)
 {
  double PI=acos(-1.0e0);
  double collision=1e-1;
  double kcon=(1e0-beam.Poisson*beam.Poisson)/beam.Young+(1e0-beam.Poisson*beam.Poisson)/1e6;
  kcon=2.5e-1*acos(-1.0e0)*1e0/kcon;
  double vcon=-2.0e0*log(collision)*sqrt(beam.Ma[ic]*kcon/(PI*PI+log(collision)*log(collision)));

  double x1[3],x2[3],normal[3];
  for(int j=0;j<3;j++) normal[j]=x[ic][j];

  double tmp=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if(tmp<1e-15){
    printf("contact wall detection error\n");
    exit(1);
  }
  for(int j=0;j<3;j++) normal[j]=normal[j]/tmp;
  double deflection=tmp+beam.rad-Radius;
  double fnormal=kcon*deflection;

  double va1[3],va2[3],vn[3];
  for(int j=0;j<3;j++) vn[j]=-v[ic][j];

  double vtmp=vn[0]*normal[0]+vn[1]*normal[1]+vn[2]*normal[2];
  double vnormal=vcon*fabs(vtmp);

  for(int j=0;j<3;j++){
    beam.fcon[ic][j]-=fnormal*normal[j];
    beam.fcon[ic][j]-=vnormal*normal[j];
  }

  //tangential force
  double mu_static=1e-1,mu_slip=9e-1*mu_static,k_static=kcon*1e-1;
  double vs=1e-3;
  double veps=1e-8*vs;
  double vt[3],ft[3],ftmp;
  //relative veloity
  for(int j=0;j<3;j++) vt[j]=vn[j]-vtmp*normal[j];
  double vt_abs=sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
  double ft_abs=mu_slip*fnormal;

  if(vt_abs>=vs){
    for(int j=0;j<3;j++) ft[j]=ft_abs*vt[j]/(vt_abs+veps);
  }else{
    for(int j=0;j<3;j++) ft[j]=k_static*vt[j]*dt+1e-1*k_static*vt[j];
    ftmp=sqrt(ft[0]*ft[0]+ft[1]*ft[1]+ft[2]*ft[2]);
      // printf("ftmp=%e,%e\n",ftmp,mu_static*fnormal);
    if(ftmp>mu_static*fnormal){
      for(int j=0;j<3;j++) ft[j]=ft_abs*vt[j]/(vt_abs+veps);
    }
  }

  for(int j=0;j<3;j++) beam.fcon[ic][j]-=ft[j];

 }

// #########################################################
/**
 * @brief contact force
 */
void multipleBeamSimulator::closedSurfaceContact(const int ibeam)
{
  double Radius=1.5e0;
  double distance,s,t;

  for(int ic=0;ic<beam[ibeam].numOfNode;ic++){
    if(beam[ibeam].mask1[ic][0]==0) continue;
    distance=sqrt(beam[ibeam].x[ic][0]*beam[ibeam].x[ic][0]
          +beam[ibeam].x[ic][1]*beam[ibeam].x[ic][1]
          +beam[ibeam].x[ic][2]*beam[ibeam].x[ic][2]);
      if(distance+beam[ibeam].rad>Radius){
        calc_wallContactForce(beam[ibeam],beam[ibeam].x,beam[ibeam].v,ic,Radius);
      }
  }
}
