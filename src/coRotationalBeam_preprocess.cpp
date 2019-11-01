#include "coRotationalBeam.h"
using namespace std;


// #################################################################
/**
 * @brief initialize physical values
 */
void coRotationalBeam::initializePhysicalValues()
{
  allocateArrays();

  for(int i=0;i<numOfElm;i++){
    l0[i] = sqrt(pow(x_ref[ie[i][0]][0]-x_ref[ie[i][1]][0],2)
           +pow(x_ref[ie[i][0]][1]-x_ref[ie[i][1]][1],2)
           +pow(x_ref[ie[i][0]][2]-x_ref[ie[i][1]][2],2));
  }

  set_mass();
  set_inertia_moment();

  set_R0();

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      u[i][j] =0e0;
      v[i][j] = 0e0;
      q[i][j] = 0e0;
      qv[i][j] = 0e0;
      udot2[i][j] = 0e0;
      qdot2[i][j] = 0e0;
    }
  }

  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<3;k++){
        ql[i][j][k] = 0e0;
        ql0[i][j][k] = 0e0;
      }
    }
  }
  set_initial_q();

  //coil parameter
  double shearRatio=G;
  //double shearRatio=Young*(2e0*(1e0+Poisson));
  // double pitch=1.1e0*wireD;
  // double n=l0[0]/pitch;
  // double D=2e0*rad-wireD;
  // EA=shearRatio*pow(wireD,4e0)*pitch/(8e0*pow(D,3e0));
  // EI=Young*pow(wireD,4e0)*pitch/(32e0*D*(2e0+Poisson));
  // GJ=Young*pow(wireD,4e0)*pitch/(64e0*D);
  // double Mass=rho*acos(-1.0e0)*pow(wireD,2e0)*acos(-1.0e0)*D*n;
  // rho=Mass/(acos(-1.0e0)*rad*rad*l0[0]);

  Poisson=(Young-2e0*shearRatio)/(2e0*shearRatio);
  EA = Young*((acos(-1.0e0)*wireD*wireD)/4e0);
  EI = Young*((acos(-1.0e0)*wireD*wireD*wireD*wireD)/64e0);
  GJ = shearRatio*((acos(-1.0e0)*wireD*wireD*wireD*wireD)/32e0);

  printf("FDstentParameter Poisson=%e EA=%e EI=%e GJ=%e rho=%e\n", Poisson,EA,EI,GJ,rho);
}

// #################################################################
/**
 * @brief allocate arrays
 */
void coRotationalBeam::allocateArrays()
{
  f=Allocation::allocate2dDOUBLE(numOfNode,3);
  T=Allocation::allocate2dDOUBLE(numOfNode,3);

  fcon=Allocation::allocate2dDOUBLE(numOfNode,3);

  f_wall=Allocation::allocate2dDOUBLE(numOfNode,3);
  fn_wall=Allocation::allocate2dDOUBLE(numOfNode,3);
  ft_wall=Allocation::allocate2dDOUBLE(numOfNode,3);
  f_contact=Allocation::allocate2dDOUBLE(numOfNode,3);

  R0=Allocation::allocate3dDOUBLE(numOfElm,3,3);

  fl=Allocation::allocate2dDOUBLE(numOfElm,7);
  fg=Allocation::allocate2dDOUBLE(numOfElm,12);

  u_local=Allocation::allocate1dDOUBLE(numOfElm);

  ql=Allocation::allocate3dDOUBLE(numOfElm,2,3);
  ql0=Allocation::allocate3dDOUBLE(numOfElm,2,3);

  ln=Allocation::allocate1dDOUBLE(numOfElm);
  l0=Allocation::allocate1dDOUBLE(numOfElm);

  u=Allocation::allocate2dDOUBLE(numOfNode,3);
  v=Allocation::allocate2dDOUBLE(numOfNode,3);
  udot2=Allocation::allocate2dDOUBLE(numOfNode,3);
  q=Allocation::allocate2dDOUBLE(numOfNode,3);
  qv=Allocation::allocate2dDOUBLE(numOfNode,3);
  qdot2=Allocation::allocate2dDOUBLE(numOfNode,3);

  energy_longitudinal=Allocation::allocate1dDOUBLE(numOfElm);
  energy_torsion=Allocation::allocate1dDOUBLE(numOfElm);
  energy_bending=Allocation::allocate1dDOUBLE(numOfElm);

  Ma=Allocation::allocate1dDOUBLE(numOfNode);
  Ia=Allocation::allocate2dDOUBLE(numOfNode,3);

  allocateArraysForAPCscheme();

  //r=Allocation::allocate1dDOUBLE(numOfNode);
  kene=Allocation::allocate1dDOUBLE(numOfNode);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      f[i][j]=0e0;
      T[i][j]=0e0;
      fcon[i][j]=0e0;
      f_wall[i][j]=0e0;
      fn_wall[i][j]=0e0;
      ft_wall[i][j]=0e0;
      f_contact[i][j]=0e0;
    }
    Ma[i]=0e0;
  }

  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<7;j++){
      fl[i][j]=0e0;
    }
    energy_longitudinal[i]=0e0;
    energy_torsion[i]=0e0;
    energy_bending[i]=0e0;
  }
}

// #################################################################
/**
 * @brief allocate arrays for APC scheme
 */
void coRotationalBeam::allocateArraysForAPCscheme()
{
  //for APC scheme only
  x_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  t_p=Allocation::allocate3dDOUBLE(numOfNode,3,3);
  u_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  q_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  v_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  qv_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  udot2_p=Allocation::allocate2dDOUBLE(numOfNode,3);
  qdot2_p=Allocation::allocate2dDOUBLE(numOfNode,3);
}
// #################################################################
/**
 * @brief initial seetings of mass and inertia tensors
 */
void coRotationalBeam::set_initial_q()
{
  double Ri[2][3][3],Rbar[2][3][3],p1[3],p2[3],Rr[3][3];

  #pragma omp parallel for private(Ri,Rbar,Rr,p1,p2)
  for(int ic=0;ic<numOfElm;ic++){

    //calc rotational matrixes
    set_Ri(Ri,t_ref,ic);
    set_Rr(Rr,x_ref,l0,Ri,p1,p2,ic);
    set_Rbar(Rbar,Rr,Ri,ic);
    //calc current displacement and rotation
    calc_theta(ql0,Rbar,ic);
  }
}
// #################################################################
/**
 * @brief set initial R0
 */
void coRotationalBeam::set_R0()
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<3;i++){
   for(int j=0;j<3;j++) R0[ic][i][j] = t_ref[ic][j][i];
    }
  }
}

// #################################################################
/**
 * @brief set mass matrrix (Vetter et al., European Journal of Mechanics A/Solids 37 (2013), 160-171)
 */
void coRotationalBeam::set_mass()
{
  for(int i=0;i<numOfElm;i++){
    Ma[ie[i][0]]  += 5e-1 * rho * PI*rad*rad * l0[i];
    Ma[ie[i][1]]  += 5e-1 * rho * PI*rad*rad * l0[i];
  }
}

// #################################################################
/**
 * @brief set initial inertia tensor (Vetter et al., European Journal of Mechanics A/Solids 37 (2013), 160-171)
 */
void coRotationalBeam::set_inertia_moment()
{
  for(int i=0;i<numOfNode;i++){
    Ia[i][0] = 4e-1*Ma[i]*pow(rad,2e0);
    Ia[i][1] = 4e-1*Ma[i]*pow(rad,2e0);
    Ia[i][2] = 4e-1*Ma[i]*pow(rad,2e0);
  }
}
