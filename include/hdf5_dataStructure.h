#ifndef _HDF5_DATA_STRUCTURE_H_
#define _HDF5_DATA_STRUCTURE_H_

//##################################################################################
//
// hdf5_dataStructure.h
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   hdf5_dataStructure.h
 * @brief  hdf5 input/output Header
 * @author T.Otani
 */

typedef struct{
  int numOfNode;   //node number
  int numOfElm;   //element number
  double rho; //density
  double mu;  //viscosity
  double rad; //radius of coil
  double wireD; //wire diamter
  double EA;  //stiffness
  double EI;  //bending stiffness
  double GJ;  //torsional stiffness
  double Young; //Young's modulus
  double Poisson; //Poisson's ratio
}exportBeamValuesHDF5;

typedef struct{
  double dt,time;
}exportValuesHDF5;

typedef struct{
  int ic1;
  int ic2;
  double ux;
  double uy;
  double uz;
}exportBeamContactHDF5;

#endif
