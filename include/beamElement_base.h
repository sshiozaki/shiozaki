//##################################################################################
//
// Multiple beam element simulator
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################

/**
 * @file   beamElement_base.h
 * @brief  beam element base class
 * @author T. Otani
 */

#ifndef _BEAMELEMENT_BASE_H_
#define _BEAMELEMENT_BASE_H_

#include "coRotationalBeam_define.h"
#include "fileIO.h"
#include "math_tools.h"
#include "TextParser.h"

class beamElement_base{
 public:
  int number;
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
  double G;     //modulus of rigidity
  double Poisson; //Poisson's ratio

  INTARRAY2 ie;      //node numbers in each element
  DOUBLEARRAY2 x,x0,x_ref;   //Position of the nodes
  DOUBLEARRAY3 t0,t_ref;    //Initial lcoal coordinates assigned at each node in each elements
  DOUBLEARRAY3 t;     //Current lcoal coordinates assigned at each node in each elements

  INTARRAY2 mask1;
  INTARRAY2 mask2;

  DOUBLEARRAY2 externalForce;    //External force of the nodes
  DOUBLEARRAY2 externalTorque;    //External torque of the nodes

  INTARRAY2 contact_beam2;
  DOUBLEARRAY1 wire_cross;

  void read_beam_geometry(TextParser &tp,const int number);

  //catheter------------------------------------------------------------------
  void read_beam_geometry_catheter(TextParser &tp, const int number);
  //--------------------------------------------------------------------------



 private:
  void read_parameters(TextParser &tp,const std::string &base);
  void setBoundaryCondition(TextParser &tp,const std::string &base);
  void set_initial_t(DOUBLEARRAY3 &t,const DOUBLEARRAY2 &x);

  void readGeometryFromFile(TextParser &tp,const std::string &base);
  void set_physicalCondition(TextParser &tp,const std::string &base);
  void tensor_initialize();

};

#endif //BEAMELMENT_BASE 
