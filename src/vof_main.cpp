//##################################################################################
//
// VOF converter
//
// Copyright (c) 2016 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################

/**
 * @file main.cpp
 * @brief VOF converter
 * @author T. Otani
 */

#include "vof_conv.h"

//return; 0 -normal
//        1 -otherwise
void VOFC::vof_main_2(const int count, const int countmax)
{

//call instance
  VOFC VOFC;

  std::string input_file;
  input_file = "minimal_test.tp"; 
  int ierror;
  if ((ierror = VOFC.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\tError at reading '%s' file\n", input_file.c_str());
    exit(1);
  }

  VOFC.initialize();

  //VOFC.make_stlfile(count, countmax);    //ファイル名設定済み　initialize()省略

  //for single mesh only
  VOFC.VOF_readGeometry_2(count);

    //calc SDF and VOF to multiDomain (all subdomain)
  int numOfSubDomain=0;

  for(unsigned int iz=0;iz<VOFC.dinfo.MPIsubDivision[2];iz++){
    for(unsigned int iy=0;iy<VOFC.dinfo.MPIsubDivision[1];iy++){
      for(unsigned int ix=0;ix<VOFC.dinfo.MPIsubDivision[0];ix++){

        VOFC.subdinfo.setSubDomain(VOFC.dinfo,ix,iy,iz);

        VOFC.VOF_converter();
        VOFC.VOF_outputData_2(numOfSubDomain, count);
        numOfSubDomain++;
      }
    }
  }

  if(numOfSubDomain>1) VOFC.exportPVTI(signedDistanceFunction);
  if(numOfSubDomain>1) VOFC.exportPVTI(volumeFraction);
  
  //calc boundary conditions

  //---------------------------post processing-----------------------------
  // std::string output_file = "VOF_" + input_file;
  // VOFC.tp.write(output_file);
  //return 0;
}

void VOFC::vof_main(){

//call instance
  VOFC VOFC;

  // if(!strcasecmp(argv[1],"--version")){
  //   std::cout << "VOF_converter version "  << VOF_VERS << std::endl;
  //   exit (0);
  // }

  // //input argument
  // if(argc!=2){
  //   printf("Invalid input. Please set tp file\n");
  //   return 1;
  // }

  // //input and check tp file
  // std::string input_file = argv[1];
  // int ierror;
  // if ((ierror = VOFC.tp.read(input_file)) != TP_NO_ERROR) {
  //  printf("\tError at reading '%s' file\n", input_file.c_str());
  //   return 1;
  // }

  std::string input_file;
  input_file = "minimal_test.tp"; 
  int ierror;
  if ((ierror = VOFC.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\tError at reading '%s' file\n", input_file.c_str());
    exit(1);
  }


  //preprocessing
  VOFC.initialize();

  //for single mesh only
  VOFC.VOF_readGeometry();

    //calc SDF and VOF to multiDomain (all subdomain)
  int numOfSubDomain=0;

  for(unsigned int iz=0;iz<VOFC.dinfo.MPIsubDivision[2];iz++){
    for(unsigned int iy=0;iy<VOFC.dinfo.MPIsubDivision[1];iy++){
      for(unsigned int ix=0;ix<VOFC.dinfo.MPIsubDivision[0];ix++){

        VOFC.subdinfo.setSubDomain(VOFC.dinfo,ix,iy,iz);

        VOFC.VOF_converter();
        VOFC.VOF_outputData(numOfSubDomain);
        numOfSubDomain++;
      }
    }
  }

  if(numOfSubDomain>1) VOFC.exportPVTI(signedDistanceFunction);
  if(numOfSubDomain>1) VOFC.exportPVTI(volumeFraction);
  //calc boundary conditions

  //---------------------------post processing-----------------------------
  std::string output_file = "VOF_" + input_file;
  VOFC.tp.write(output_file);
  //return 0;

}
