/**
 * @file coRotationalBeamElasticForce.cpp
 * @brief coRotationalBeam class
 * @author T. Otani
 */

#include "coRotationalBeam.h"
using namespace std;

// #################################################################
/**
 * @brief calc elastic force
 * @param [in] x_p    node position vector (predictor)
 * @param [in] t_p    unit basis (predictor)
 */
void coRotationalBeam::calc_elastic_force(const DOUBLEARRAY2 &x_p,const DOUBLEARRAY3 &t_p)
{
  double Ri[2][3][3],Rbar[2][3][3],Rr[3][3],p1[3],p2[3],B[7][12],F[12];

  // printf("localVaribales\n");
  calc_current_length(x_p);

  #pragma omp parallel for private(Ri,Rbar,Rr,p1,p2,B,F)
  for(int ic=0;ic<numOfElm;ic++){

    //calc rotational matrixes
    set_Ri(Ri,t_p,ic);
    set_Rr(Rr,x_p,ln,Ri,p1,p2,ic);
    set_Rbar(Rbar,Rr,Ri,ic);

    //calc current displacement and rotation
    calc_theta(ql,Rbar,ic);
    u_local[ic] = ln[ic]-l0[ic];

    for(int p=0;p<2;p++){
      for(int j=0;j<3;j++){
        ql[ic][p][j]-=ql0[ic][p][j];
      }
    }
    //calc force in each lcoal coordinates
    calc_force_in_local_coordinates(ic);

    //calc elastic energy (longitudinal, torsion and bending)
    calc_elastic_energy(ic);

    //calc translation matrix
    calc_BMatrix(B,Rr,p1,p2,ic);

    //calc force in glboal coordinates
    calc_force_in_global_coordinates(F,B,ic);

    //transform from spatial rotation to rotation vector
    set_incrementalRotationalVector(ic,F);
  }
  calc_nodal_force();
}

// #################################################################
/**
 * @brief calc elastic force
 */
void coRotationalBeam::calc_elastic_force()
{
  double Ri[2][3][3],Rbar[2][3][3],Rr[3][3],p1[3],p2[3],B[7][12],F[12];

  // printf("localVaribales\n");
  calc_current_length(x);

  #pragma omp parallel for private(Ri,Rbar,Rr,p1,p2,B,F)
  for(int ic=0;ic<numOfElm;ic++){

    //calc rotational matrixes
    set_Ri(Ri,t,ic);
    set_Rr(Rr,x,ln,Ri,p1,p2,ic);
    set_Rbar(Rbar,Rr,Ri,ic);

    //calc current displacement and rotation
    calc_theta(ql,Rbar,ic);
    u_local[ic] = ln[ic]-l0[ic];
    for(int p=0;p<2;p++){
      for(int j=0;j<3;j++){
        ql[ic][p][j]-=ql0[ic][p][j];
      }
    }
    //calc force in each lcoal coordinates
    calc_force_in_local_coordinates(ic);

    //calc translation matrix
    calc_BMatrix(B,Rr,p1,p2,ic);

    //calc force in glboal coordinates
    calc_force_in_global_coordinates(F,B,ic);

    //transform from spatial rotation to rotation vector
    set_incrementalRotationalVector(ic,F);
  }
  calc_nodal_force();
}

// #################################################################
/**
 * @brief calc nodal force
 */
void coRotationalBeam::calc_nodal_force()
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0;i<3;i++){
      f[ic][i] = 0e0;
      T[ic][i] = 0e0;
    }
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<3;i++){
      f[ie[ic][0]][i] += fg[ic][i];
      T[ie[ic][0]][i] += fg[ic][i+3];
      f[ie[ic][1]][i] += fg[ic][i+6];
      T[ie[ic][1]][i] += fg[ic][i+9];
    }
  }
}

// #################################################################
/**
 * @brief set incremental rotation vector
 * @param [in] ic     element number
 * @param [in] F      force vector in global coordinates
 */
void coRotationalBeam::set_incrementalRotationalVector(const int &ic,const double (&F)[12])
{
  double Br[12][12],tmp,T[3][3];
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      Br[i][j]=0e0;
    }
    Br[i][i]=1e0;
  }

  tmp=sqrt(pow(q[ie[ic][0]][0],2e0)+pow(q[ie[ic][0]][1],2e0)+pow(q[ie[ic][0]][2],2e0));
  if(tmp>1e-20){
    calc_T(T,ie[ic][0],tmp);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Br[i+3][j+3]=T[i][j];
      }
    }
  }

  tmp=sqrt(pow(q[ie[ic][1]][0],2e0)+pow(q[ie[ic][1]][1],2e0)+pow(q[ie[ic][1]][2],2e0));
  if(tmp>1e-20){
    calc_T(T,ie[ic][1],tmp);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Br[i+9][j+9]=T[i][j];
      }
    }
  }

  for(int i=0;i<12;i++){
    fg[ic][i]=0e0;
    for(int j=0;j<12;j++){
      fg[ic][i]+=Br[j][i]*F[j];
    }
  }

}

// #################################################################
/**
 * @brief calc T
 * @param [out]    T  Temporal coordinate transformation matrix
 * @param [in]    ic  element number
 * @param [in] theta  rotation angle
 */
void coRotationalBeam::calc_T(double (&T)[3][3],const int &ic,const double &theta)
{
  double I[3][3],S[3][3],S2[3][3];
  I[0][0]=1e0; I[0][1]=0e0; I[0][2]=0e0;
  I[1][0]=0e0; I[1][1]=1e0; I[1][2]=0e0;
  I[2][0]=0e0; I[2][1]=0e0; I[2][2]=1e0;

  S[0][0]=0e0; S[0][1]=-q[ic][2]; S[0][2]=q[ic][1];
  S[1][0]=q[ic][2]; S[1][1]=0e0; S[1][2]=-q[ic][0];
  S[2][0]=-q[ic][1]; S[2][1]=q[ic][0]; S[2][2]=0e0;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      S2[i][j]=0e0;
      for(int k=0;k<3;k++){
        S2[i][j]+=S[i][k]*S[k][j];
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      T[i][j]=I[i][j]+(1e0-cos(theta))/(theta*theta)*S[i][j]+(theta-sin(theta))/(theta*theta*theta)*S2[i][j];
    }
  }

}

// #################################################################
/**
 * @brief calc elastic force (glboal coordinates)
 * @param [out] F     force vector in global coordinates
 * @param [in]  B     Coordinate transformation matrix
 * @param [in]  ic    element number
 */
void coRotationalBeam::calc_force_in_global_coordinates(double (&F)[12],const double (&B)[7][12],const int &ic)
{
  for(int i=0;i<12;i++){
    F[i] = 0e0;
    for(int j=0;j<7;j++) F[i] += B[j][i] * fl[ic][j];
  }
}
// #################################################################
/**
 * @brief calc B matrix
 * @param [out]  B      Coordinate transformation matrix
 * @param [out]  Rr     Rotation matrix [global to current (element level)]
 * @param [in]  p1      indicator vector
 * @param [in]  p2      indicator vector
 * @param [in]  ic      element number
 */
void coRotationalBeam::calc_BMatrix(double (&B)[7][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic)
{
  double Ba[7][7],Bg[7][12];

  calc_Ba(Ba,ic);
  calc_Bg(Bg,Rr,p1,p2,ic);

  for(int i=0;i<7;i++){
    for(int j=0;j<12;j++){
      B[i][j] = 0e0;
      for(int k=0;k<7;k++) B[i][j] += Ba[i][k] * Bg[k][j];
    }
  }
}
// #################################################################
/**
 * @brief calc Bg
 * @param [out]  Bg     Coordinate transformation matrix
 * @param [out]  Rr     Rotation matrix [global to current (element level)]
 * @param [in]   p1     indicator vector
 * @param [in]   p2     indicator vector
 * @param [in]   ic     element number
 */
void coRotationalBeam::calc_Bg(double (&Bg)[7][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic)
{
  double P[6][12],E[12][12],Id[6][12],trG[3][12];

  for(int i=0;i<6;i++){
    for(int j=0;j<12;j++) Id[i][j] = 0e0;
  }

  Id[0][3] = 1e0; Id[1][4] = 1e0; Id[2][5] = 1e0;
  Id[3][9] = 1e0; Id[4][10] = 1e0; Id[5][11] = 1e0;

  for(int i=0;i<7;i++){
    for(int j=0;j<12;j++) Bg[i][j] = 0e0;
  }

  Bg[0][0] = -Rr[0][0];
  Bg[0][1] = -Rr[1][0];
  Bg[0][2] = -Rr[2][0];
  Bg[0][6] = Rr[0][0];
  Bg[0][7] = Rr[1][0];
  Bg[0][8] = Rr[2][0];

  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++) E[i][j] = 0e0;
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      E[i][j]   = Rr[i][j];
      E[i+3][j+3] = Rr[i][j];
      E[i+6][j+6] = Rr[i][j];
      E[i+9][j+9] = Rr[i][j];
    }
  }

  calc_trG(trG,Rr,p1,p2,ic);

  for(int i=0;i<3;i++){
    for(int j=0;j<12;j++){
      P[i][j]  = Id[i][j]  - trG[i][j];
      P[i+3][j] = Id[i+3][j] - trG[i][j];
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<12;j++){
      for(int k=0;k<12;k++) Bg[i+1][j] += P[i][k] * E[j][k];
    }
  }
}
// #################################################################
/**
 * @brief calc trG
 * @param [out] trG   transpose of G matrix
 * @param [in]  Rr    Rotation matrix [global to current (element level)]
 * @param [in] p1     indicator vector
 * @param [in] p2     indicator vector
 * @param [in] ic     element number
 */
void coRotationalBeam::calc_trG(double (&trG)[3][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic)
{
  double p_1,p_2,p_11,p_12,p_13,p_21,p_22,p_23,e,e11,e12,e21,e22;
  double p[3];

  for(int i=0;i<3;i++) p[i] = 5e-1 * (p1[i] + p2[i]);

  p_1 = Rr[0][0] * p[0] + Rr[1][0] * p[1] + Rr[2][0] * p[2];
  p_2 = Rr[0][1] * p[0] + Rr[1][1] * p[1] + Rr[2][1] * p[2];

  p_11 = Rr[0][0] * p1[0] + Rr[1][0] * p1[1] + Rr[2][0] * p1[2];
  p_12 = Rr[0][1] * p1[0] + Rr[1][1] * p1[1] + Rr[2][1] * p1[2];
  p_13 = Rr[0][2] * p1[0] + Rr[1][2] * p1[1] + Rr[2][2] * p1[2];

  p_21 = Rr[0][0] * p2[0] + Rr[1][0] * p2[1] + Rr[2][0] * p2[2];
  p_22 = Rr[0][1] * p2[0] + Rr[1][1] * p2[1] + Rr[2][1] * p2[2];
  p_23 = Rr[0][2] * p2[0] + Rr[1][2] * p2[1] + Rr[2][2] * p2[2];

  e   = p_1 /p_2;
  e11 = p_11/p_2;
  e12 = p_12/p_2;
  e21 = p_21/p_2;
  e22 = p_22/p_2;

  for(int i=0;i<3;i++){
    for(int j=0;j<12;j++) trG[i][j] = 0e0;
  }

  trG[0][2]  =  e / ln[ic];
  trG[0][3]  =  5e-1*e12;
  trG[0][4]  = -5e-1*e11;
  trG[0][8]  = -e / ln[ic];
  trG[0][9]  =  5e-1*e22;
  trG[0][10] = -5e-1*e21;
  trG[1][2]  =  1e0 / ln[ic];
  trG[1][8]  = -1e0 / ln[ic];
  trG[2][1]  = -1e0 / ln[ic];
  trG[2][7]  =  1e0 / ln[ic];
}
// #################################################################
/**
 * @brief calc Ba
 * @param [out] Ba  Coordinate transformation matrix
 * @param [in] ic   element number
 */
void coRotationalBeam::calc_Ba(double (&Ba)[7][7],const int &ic)
{
  double q0,tmp,invT[2][3][3];

  for(int i=0;i<7;i++){
    for(int j=0;j<7;j++) Ba[i][j] = 0e0;
  }

  for(int p=0;p<2;p++){

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) invT[p][i][j] = 0e0;
    }

    q0=sqrt(pow(ql[ic][p][0],2e0)+pow(ql[ic][p][1],2e0)+pow(ql[ic][p][2],2e0));

    if(q0>1e-20){
      tmp=5e-1*q0/tan(5e-1*q0);
      for(int i=0;i<3;i++) invT[p][i][i] += tmp;

      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          invT[p][i][j] += (1e0-tmp)/(q0*q0) * ql[ic][p][i]*ql[ic][p][j];
        }
      }

      invT[p][0][1] -= -5e-1 * ql[ic][p][2];
      invT[p][0][2] -= 5e-1 * ql[ic][p][1];
      invT[p][1][0] -= 5e-1 * ql[ic][p][2];
      invT[p][1][2] -= -5e-1 * ql[ic][p][0];
      invT[p][2][0] -= -5e-1 * ql[ic][p][1];
      invT[p][2][1] -= 5e-1 * ql[ic][p][0];

    }else{
      invT[p][0][0] = 1e0;
      invT[p][1][1] = 1e0;
      invT[p][2][2] = 1e0;
    }
  }

  Ba[0][0] = 1e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Ba[i+1][j+1] = invT[0][i][j];
      Ba[i+4][j+4] = invT[1][i][j];
    }
  }
}

// #################################################################
/**
 * @brief calc elastic energy
 */
void coRotationalBeam::calc_elastic_energy(const int &ic)
{
  energy_longitudinal[ic] = 5e-1*EA/l0[ic]*u_local[ic]*u_local[ic];
  energy_torsion[ic] = 5e-1*GJ/l0[ic]*(ql[ic][0][0]-ql[ic][1][0])*(ql[ic][0][0]-ql[ic][1][0]);
  energy_bending[ic] = 2e0*EI/l0[ic]*(ql[ic][0][1]*ql[ic][0][1]+ql[ic][0][1]*ql[ic][1][1]+ql[ic][1][1]*ql[ic][1][1]
                                      +ql[ic][0][2]*ql[ic][0][2]+ql[ic][0][2]*ql[ic][1][2]+ql[ic][1][2]*ql[ic][1][2]);
}


// #################################################################
/**
 * @brief calc elastic force (local coordinates)
 * @param [in] ic element number
 */
void coRotationalBeam::calc_force_in_local_coordinates(const int &ic)
{
  fl[ic][0] = EA / l0[ic] * u_local[ic];
  fl[ic][1] = GJ / l0[ic] * (ql[ic][0][0]-ql[ic][1][0]);
  fl[ic][2] = EI / l0[ic] * ( 4e0*ql[ic][0][1] + 2e0 * ql[ic][1][1] );
  fl[ic][3] = EI / l0[ic] * ( 4e0*ql[ic][0][2] + 2e0 * ql[ic][1][2] );
  fl[ic][4] = GJ / l0[ic] * (  -ql[ic][0][0] + ql[ic][1][0] );
  fl[ic][5] = EI / l0[ic] * ( 2e0*ql[ic][0][1] + 4e0 * ql[ic][1][1] );
  fl[ic][6] = EI / l0[ic] * ( 2e0*ql[ic][0][2] + 4e0 * ql[ic][1][2] );

  // fl[ic][1] = GJ / l0[ic] * ( (ql[ic][0][0]-ql0[ic][0][0]) - (ql[ic][1][0]-ql0[ic][1][0]) );
  // fl[ic][2] = EI / l0[ic] * ( 4e0*(ql[ic][0][1]-ql0[ic][0][1]) + 2e0 * (ql[ic][1][1]-ql0[ic][1][1]) );
  // fl[ic][3] = EI / l0[ic] * ( 4e0*(ql[ic][0][2]-ql0[ic][0][2]) + 2e0 * (ql[ic][1][2]-ql0[ic][1][2]) );
  // fl[ic][4] = GJ / l0[ic] * (  -(ql[ic][0][0]-ql0[ic][0][0]) + (ql[ic][1][0]-ql0[ic][1][0]) );
  // fl[ic][5] = EI / l0[ic] * ( 2e0*(ql[ic][0][1]-ql0[ic][0][1]) + 4e0 * (ql[ic][1][1]-ql0[ic][1][1]) );
  // fl[ic][6] = EI / l0[ic] * ( 2e0*(ql[ic][0][2]-ql0[ic][0][2]) + 4e0 * (ql[ic][1][2]-ql0[ic][1][2]) );
}
// #################################################################
/**
 * @brief calc theta (Nour-Omid and Rankin, Compt. Methods Appl. Mech. Eng., 1991.)
 * @param [out] ql      rotation angle in local coordinates
 * @param [in] Rbar     rotation matrix [reference to current coorindate (node level)]
 * @param [in] ic element number
 */
void coRotationalBeam::calc_theta(DOUBLEARRAY3 &ql,const double (&Rbar)[2][3][3],const int &ic)
{
  double trR,Ra[3][3],tau;

  for(int p=0;p<2;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Ra[i][j] = Rbar[p][i][j]-Rbar[p][j][i];
      }
    }
    ql[ic][p][0] = Ra[2][1];
    ql[ic][p][1] = Ra[0][2];
    ql[ic][p][2] = Ra[1][0];
    tau=5e-1*sqrt(ql[ic][p][0]*ql[ic][p][0]+ql[ic][p][1]*ql[ic][p][1]+ql[ic][p][2]*ql[ic][p][2]);
    if(tau<1e-20){
      for(int i=0;i<3;i++) ql[ic][p][i]=0e0;
    }else{
      for(int i=0;i<3;i++) ql[ic][p][i]=5e-1*asin(tau)/tau*ql[ic][p][i];
    }
    // for(int i=0;i<3;i++){
    //   for(int j=0;j<3;j++)Ra[i][j]*=5e-1;
    // }
    // trR=Rbar[p][0][0]+Rbar[p][1][1]+Rbar[p][2][2];
    // ql[ic][p][0] = 4e0*Ra[2][1]/(1e0+trR);
    // ql[ic][p][1] = 4e0*Ra[0][2]/(1e0+trR);
    // ql[ic][p][2] = 4e0*Ra[1][0]/(1e0+trR);

  }
}
// #################################################################
/**
 * @brief calc Rbar
 * @param [out] Rbar     rotation matrix [reference to current coorindate (node level)]
 * @param [in] Rr       rotation matrix [global to current coorindate (element level)]
 * @param [in] Ri       rotation matrix [global to reference coorindate (node level)]
 * @param [in] ic       element number
 */
void coRotationalBeam::set_Rbar(double (&Rbar)[2][3][3],const double (&Rr)[3][3],const double (&Ri)[2][3][3],const int &ic)
{
  for(int p=0;p<2;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Rbar[p][i][j]=0e0;
        for(int k=0;k<3;k++){
          for(int l=0;l<3;l++) Rbar[p][i][j] += Rr[k][i] * Ri[p][k][l] * R0[ic][l][j]; 
        }
      }
    }
  }
}
// #################################################################
/**
 * @brief calc Rr
 * @param [out] Rr        rotation matrix [global to current coorindate (element level)]
 * @param [in] ln         element length in current configuration
 * @param [in] Ri         rotation matrix [global to reference coorindate (node level)]
 * @param [out] p1        indicator vector
 * @param [out] p2        indicator vector
 * @param [ic] ic         element number
 */
void coRotationalBeam::set_Rr(double (&Rr)[3][3],const DOUBLEARRAY2 &x,const DOUBLEARRAY1 &ln,const double (&Ri)[2][3][3],double (&p1)[3],double (&p2)[3],const int &ic)
{
  double r1[3],r2[3],r3[3],p[3],dr2,dr3;

  //calc r1
  for(int i=0;i<3;i++) r1[i] = (x[ie[ic][1]][i]-x[ie[ic][0]][i]) / ln[ic];
  double tmp=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
  for(int i=0;i<3;i++) r1[i]/=tmp;

  //calc p
  for(int i=0;i<3;i++){
    p1[i] = 0e0;
    p2[i] = 0e0;
    for(int j=0;j<3;j++){
      p1[i] += Ri[0][i][j] * R0[ic][j][1];
      p2[i] += Ri[1][i][j] * R0[ic][j][1];
    }
  }
  for(int i=0;i<3;i++) p[i]=5e-1*(p1[i]+p2[i]);

  //calc r2 and r3
  mathTool::crossProduct(r1,p,r3,dr3);
  for(int i=0;i<3;i++) r3[i]=r3[i]/dr3;
  mathTool::crossProduct(r3,r1,r2,dr2);

  for(int i=0;i<3;i++){
    Rr[i][0] = r1[i];
    Rr[i][1] = r2[i];
    Rr[i][2] = r3[i];
  }
}
// #################################################################
/**
 * @brief calc Ri
 * @param [out] Ri         rotation matrix [global to reference coorindate (node level)]
 * @param [in]  t          unit basis
 * @param [in] ic          element number
 */
void coRotationalBeam::set_Ri(double (&Ri)[2][3][3],const DOUBLEARRAY3 &t,const int &ic)
{
  for(int p=0;p<2;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Ri[p][i][j] = 0e0;
        for(int k=0;k<3;k++){
          // Ri[p][i][j] += t[ie[ic][p]][k][i] * t0[ic][k][j];
          Ri[p][i][j] += t[ie[ic][p]][k][i] * t_ref[ic][k][j];
        }
      }
    }
  }
}
