#ifndef _COROTATIONALBEAM_DEFINE_H_
#define _COROTATIONALBEAM_DEFINE_H_

//##################################################################################
//
// COROTATIONALBEAM Base
//
// Copyright (c) 2017 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file  coRotationalBeam_define.h
 * @brief  Definition Header
 * @author T. Otani
 */
#include <string>
#include <vector>
// #include <boost/multi_array.hpp>
#include "allocation.h"

#define GRAVITY 9.80665 // (m/s2)
static const double PI = 3.1415926535897932384626;

// typedef std::vector<int> INTVECTOR1;
// typedef std::vector<std::vector<int>> INTVECTOR2;
// typedef std::vector<std::vector<std::vector<int>>> INTVECTOR3;
// typedef std::vector<std::vector<std::vector<std::vector<int>>>> INTVECTOR4;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> INTVECTOR5;

// typedef std::vector<double> DOUBLEVECTOR1;
// typedef std::vector<std::vector<double>> DOUBLEVECTOR2;
// typedef std::vector<std::vector<std::vector<double>>> DOUBLEVECTOR3;
// typedef std::vector<std::vector<std::vector<std::vector<double>>>> DOUBLEVECTOR4;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> DOUBLEVECTOR5;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DOUBLEVECTOR6;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> DOUBLEVECTOR7;

// typedef boost::multi_array<int,1> INTARRAY1;
// typedef boost::multi_array<int,2> INTARRAY2;
// typedef boost::multi_array<int,3> INTARRAY3;
// typedef boost::multi_array<int,4> INTARRAY4;
// typedef boost::multi_array<int,5> INTARRAY5;

// typedef boost::multi_array<double,1> DOUBLEARRAY1;
// typedef boost::multi_array<double,2> DOUBLEARRAY2;
// typedef boost::multi_array<double,3> DOUBLEARRAY3;
// typedef boost::multi_array<double,4> DOUBLEARRAY4;
// typedef boost::multi_array<double,5> DOUBLEARRAY5;
// #else
//   typedef int*        INTARRAY1;
//   typedef int**       INTARRAY2;
//   typedef int***      INTARRAY3;
//   typedef int****     INTARRAY4;
//   typedef int*****    INTARRAY5;
//   typedef double*     DOUBLEARRAY1;
//   typedef double**    DOUBLEARRAY2;
//   typedef double***   DOUBLEARRAY3;
//   typedef double****  DOUBLEARRAY4;
//   typedef double***** DOUBLEARRAY5;
// #endif

// general
#define ON      1
#define OFF      0
#define GMM      1
#define observer    2

//data encode
#define INT       0
#define DOUBLE      1
#define ASCII      0
#define BINARY      1


#endif // _COROTATIONALBEAM_DEFINE_H_
