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
 * @file contactBox.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"
#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

// #########################################################
/**
 * @brief register node position in box
 */
void multipleBeamSimulator::register_box()
{
  int ix[3];
  double x[3];

  #pragma omp parallel for private(x,ix)
  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int ic=0;ic<beam[ibeam].numOfElm;ic++){

      for(int j=0;j<3;j++){
        x[j] = 5e-1*(beam[ibeam].x[beam[ibeam].ie[ic][0]][j]+beam[ibeam].x[beam[ibeam].ie[ic][1]][j]);
        ix[j] = (int)((x[j]-origin[j])/dx[j]);
        if(ix[j]<0 || numOfVoxel[j]<=ix[j]) ix[j]=-1;
        box[ibeam][ic][j] = ix[j];
      }

    }
  }

  int num;
  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int i=0;i<beam[ibeam].numOfElm;i++){
      if(box[ibeam][i][0]==-1) continue;
      if(box[ibeam][i][1]==-1) continue;
      if(box[ibeam][i][2]==-1) continue;
      for(int j=0;j<3;j++) ix[j]=box[ibeam][i][j];
      num = numOfNodeInBox[ix[0]][ix[1]][ix[2]][ibeam];
      nodeInBox[ix[0]][ix[1]][ix[2]][ibeam][num] = i;
      numOfNodeInBox[ix[0]][ix[1]][ix[2]][ibeam] += 1;
    }
  }

}

// #########################################################
/**
 * @brief delete the box information
 */
void multipleBeamSimulator::free_box()
{
  // #pragma omp parallel for
  // for(int ix=0,nx=numOfVoxel[0];ix<nx;++ix){
  //   for(int iy=0,ny=numOfVoxel[1];iy<ny;++iy){
  //     for(int iz=0,nz=numOfVoxel[2];iz<nz;++iz){
  //       for(int ibeam=0;ibeam<numOfBeams;++ibeam){
  //         for(int i=0,n=numOfNodeInBox[ix][iy][iz][ibeam];i<n;i++) nodeInBox[ix][iy][iz][ibeam][i] = 0;
  //         numOfNodeInBox[ix][iy][iz][ibeam] = 0;
  //       }
  //     }
  //   }
  // }

  #pragma omp parallel for
  for(int ix=0;ix<numOfVoxel[0];++ix){
    for(int iy=0;iy<numOfVoxel[1];++iy){
      for(int iz=0;iz<numOfVoxel[2];++iz){
        for(int ibeam=0;ibeam<numOfBeams;++ibeam){
          for(int i=0,n=numOfNodeInBox[ix][iy][iz][ibeam];i<n;i++) nodeInBox[ix][iy][iz][ibeam][i] = 0;
          numOfNodeInBox[ix][iy][iz][ibeam] = 0;
        }
      }
    }
  }
}

// #########################################################
/**
 * @brief initialize the box
 */
void multipleBeamSimulator::initializeBox()
{
  int maxElm=0;
  double maxl0=0e0;

  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int ic=0;ic<beam[ibeam].numOfElm;ic++){
     if(maxl0<beam[ibeam].l0[ic]) maxl0=beam[ibeam].l0[ic];
    }
  }

  for(int i=0;i<3;i++){
    origin[i]=wall[0].origin[i];
    numOfVoxel[i] = (int)(wall[0].glboalLength[i]/(1.1e0*maxl0));
    dx[i] = wall[0].glboalLength[i] / (double)numOfVoxel[i];
  }
  printf("wall numberfVoxel=%d(x) %d(y) %d(z)\n",numOfVoxel[0],numOfVoxel[1],numOfVoxel[2]);

  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    if(maxElm<beam[ibeam].numOfElm) maxElm = beam[ibeam].numOfElm;
  }

  box=Allocation::allocate3dINT(numOfBeams,maxElm,3);
  numOfNodeInBox=Allocation::allocate4dINT(numOfVoxel[0],numOfVoxel[1],numOfVoxel[2],numOfBeams);
  nodeInBox=Allocation::allocate5dINT(numOfVoxel[0],numOfVoxel[1],numOfVoxel[2],numOfBeams,50);

  for(int ix=0;ix<numOfVoxel[0];ix++){
    for(int iy=0;iy<numOfVoxel[1];iy++){
      for(int iz=0;iz<numOfVoxel[2];iz++){
        for(int ibeam=0;ibeam<numOfBeams;ibeam++){
          numOfNodeInBox[ix][iy][iz][ibeam] = 0;
        }
      }
    }
  }

}


// #########################################################
/**
 * @brief initialize the box
 */
void multipleBeamSimulator::initializeBox_catheter()
{
  int maxElm=0;
  double maxl0=0e0;

  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int ic=0;ic<beam[ibeam].numOfElm;ic++){
     if(maxl0<beam[ibeam].l0[ic]) maxl0=beam[ibeam].l0[ic];
    }
  }

  for(int i=0;i<3;i++){
    origin[i]=wall[1].origin[i];
    numOfVoxel[i] = (int)(wall[1].glboalLength[i]/(1.1e0*maxl0));
    dx[i] = wall[1].glboalLength[i] / (double)numOfVoxel[i];
  }
  printf("wall numberfVoxel=%d(x) %d(y) %d(z)\n",numOfVoxel[0],numOfVoxel[1],numOfVoxel[2]);

  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    if(maxElm<beam[ibeam].numOfElm) maxElm = beam[ibeam].numOfElm;
  }

  box=Allocation::allocate3dINT(numOfBeams,maxElm,3);
  numOfNodeInBox=Allocation::allocate4dINT(numOfVoxel[0],numOfVoxel[1],numOfVoxel[2],numOfBeams);
  nodeInBox=Allocation::allocate5dINT(numOfVoxel[0],numOfVoxel[1],numOfVoxel[2],numOfBeams,50);

  for(int ix=0;ix<numOfVoxel[0];ix++){
    for(int iy=0;iy<numOfVoxel[1];iy++){
      for(int iz=0;iz<numOfVoxel[2];iz++){
        for(int ibeam=0;ibeam<numOfBeams;ibeam++){
          numOfNodeInBox[ix][iy][iz][ibeam] = 0;
        }
      }
    }
  }

}
