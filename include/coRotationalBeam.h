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

#ifndef _COROTATIONALBEAM_H_
#define _COROTATIONALBEAM_H_

/**
 * @file   coRotationalBeam.h
 * @brief  Single beam treatment main Header
 * @author T. Otani
 */

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <strings.h>
#include <sys/stat.h>
#include <omp.h>
#include <chrono>

#include "beamElement_base.h"
#include "coRotationalBeam_define.h"
#include "math_tools.h"
#include "TextParser.h"

class coRotationalBeam : public beamElement_base{
 public:

  double uref,etaMax,etaMin,etaBar;

  DOUBLEARRAY2 f;     //Force fo the nodes
  DOUBLEARRAY2 T;     //Torque fo the nodes

  //DOUBLEARRAY2 f_add;     //Force fo the nodes
  //DOUBLEARRAY2 T_add;     //Torque fo the nodes

  //DOUBLEARRAY2 spforce;
  
  DOUBLEARRAY1 u_local;  //current displacement(local coordinates)
  DOUBLEARRAY3 ql;    //current angle in element coordinates
  DOUBLEARRAY3 ql0;   //initial angle in element coordinates

  DOUBLEARRAY2 fl;    //cuurent inner force in the local coordinates (non-addmissive)
  DOUBLEARRAY2 fg;    //current inner force in the glbal coordinates
  DOUBLEARRAY2 v;     //Velocity of the nodes
  DOUBLEARRAY2 qv;    //Angular velocity of the nodes
  DOUBLEARRAY2 q;     //Rotation angle of nodes

  DOUBLEARRAY2 f_wall;
  DOUBLEARRAY2 fn_wall;
  DOUBLEARRAY2 ft_wall;
  DOUBLEARRAY2 f_contact;

  DOUBLEARRAY1 energy_longitudinal;
  DOUBLEARRAY1 energy_torsion;
  DOUBLEARRAY1 energy_bending;

  DOUBLEARRAY1 Ma;    //Mass

  DOUBLEARRAY2 fcon;
  DOUBLEARRAY2 x_p,u_p,q_p,v_p,qv_p,udot2_p,qdot2_p;  //for explicitAPC scheme
  DOUBLEARRAY3 t_p;

  DOUBLEARRAY1 l0;    //initial length
  DOUBLEARRAY1 ln;    //current length
  DOUBLEARRAY2 u;     //displacement
  DOUBLEARRAY2 udot2,qdot2;
  DOUBLEARRAY3 R0;    //Element rotation tensor from global to reference local coordinates
  DOUBLEARRAY2 Ia;    //Inertia moment (global coordinate)

  //DOUBLEARRAY1 r;
  DOUBLEARRAY1 kene;

  //DOUBLEARRAY1 ini_distance;

 private:

 //coRotationalBeam_preprocess.cpp
 public:
  void initializePhysicalValues();

 private:
  void set_mass();
  void set_inertia_moment();
  void set_initial_q();
  void set_R0();
  void allocateArrays();
  void allocateArraysForAPCscheme();

  //coRotationalBeam.cpp
 public:
  void explicitScheme();
  void explicitAPCScheme();
  void quasiStaticScheme();
  void exportVTP_polyLine(const std::string &file);
  void exportVTP_polyLine(const std::string &file,const DOUBLEARRAY2 &x1,const DOUBLEARRAY3 &t1);

  void quasiStaticScheme(const double dt);
  void semiImplicit_Euler_scheme(const double dt);
  void predictor(const double dt);
  void corrector(const double dt);
  double estimator(const double uref,const double qref,const double dt);
  void calc_udot2(const double dt);
  void calc_elastic_force();
  void calc_elastic_force(const DOUBLEARRAY2 &x_p,const DOUBLEARRAY3 &t_p);
  //void calc_t(DOUBLEARRAY3 &t,const DOUBLEARRAY2 &q1);
  void calc_t(DOUBLEARRAY3 &t,const DOUBLEARRAY3 &t_pre,const DOUBLEARRAY2 &q1);

  void calc_elastic_energy(const int &ic);

  void calc_r(DOUBLEARRAY1 &r, const DOUBLEARRAY2 &x);
  void calc_kinetic_energy(DOUBLEARRAY1 &kene, const DOUBLEARRAY2 &v);
  double calc_total(DOUBLEARRAY1 &total1);
  double calc_max(DOUBLEARRAY1 &max1, const double max2);
  

 private:
  void calc_current_inertia_moment();
  void calc_R(double (&R)[3][3],const DOUBLEARRAY2 &q1,const int &ic);
  void calc_current_length(const DOUBLEARRAY2 &x);

  void set_incrementalRotationalVector(const int &ic,const double (&F)[12]);
  void calc_T(double (&T)[3][3],const int &ic,const double &theta);


  //coRotationalBeamElasticForce.cpp
 private:
  void set_Ri(double (&Ri)[2][3][3],const DOUBLEARRAY3 &t,const int &ic);
  void set_Rr(double (&Rr)[3][3],const DOUBLEARRAY2 &x,const DOUBLEARRAY1 &ln,const double (&Ri)[2][3][3],double (&p1)[3],double (&p2)[3],const int &ic);
  void set_Rbar(double (&Rbar)[2][3][3],const double (&Rr)[3][3],const double (&Ri)[2][3][3],const int &ic);
  void calc_theta(DOUBLEARRAY3 &ql,const double (&Rbar)[2][3][3],const int &ic);
  void calc_Ba(double (&Ba)[7][7],const int &ic);
  void calc_Bg(double (&Bg)[7][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic);
  void calc_trG(double (&trG)[3][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic);
  void calc_force_in_local_coordinates(const int &ic);
  void calc_BMatrix(double (&B)[7][12],const double (&Rr)[3][3],const double (&p1)[3],const double (&p2)[3],const int &ic);
  void calc_force_in_global_coordinates(double (&F)[12],const double (&B)[7][12],const int &ic);
  void calc_nodal_force();
};

#endif //COROTATIONALBEAM
