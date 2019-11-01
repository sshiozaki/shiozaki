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
 * @file input_hdf5.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <H5Cpp.h>

using namespace std;
using namespace H5;


int multipleBeamSimulator::input_hdf5_all(const std::string &fileName,const std::string &groupName)
{
   // Try block to detect exceptions raised by any of the calls inside it
  try
  {
    Exception::dontPrint();

    H5std_string FILE_NAME(fileName.c_str());
    H5File file(FILE_NAME,H5F_ACC_RDONLY);
    string dataName;
    dataName=groupName+"/contact/u_stat_wall";
    HDF5input::read_3dDOUBLEARRAY(file,dataName,u_stat_wall);
    dataName=groupName+"/contact/wallContact";
    HDF5input::read_2dINTARRAY(file,dataName,wallContact);

    for(int ibeam=0;ibeam<numOfBeams;ibeam++){

      //read scalars

      //read arrays
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ie";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].ie);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/R0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].R0);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ql0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql0);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/l0";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].l0);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/Ma";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].Ma);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/Ia";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].Ia);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/x0";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x0);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/x_ref";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x_ref);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/t0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t0);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/t_ref";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t_ref);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/externalForce";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalForce);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/externalTorque";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalTorque);

      dataName=groupName+"/beam_"+to_string(ibeam)+"/x";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/t";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/mask1";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].mask1);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/mask2";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].mask2);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fl";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fl);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fg";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fg);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/force";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].f);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/Torque";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].T);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fcon";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fcon);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/u_local";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].u_local);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ql";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ln";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].ln);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/u";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].u);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/v";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].v);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/udot2";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].udot2);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/q";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].q);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/qv";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].qv);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/qdot2";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].qdot2);
    }
    // for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    //   groupName2=groupName+"/beam_" +to_string(ibeam);
    //   file.createGroup(groupName2.c_str());
    // }

    // //scalar value export
    // exportValuesHDF5 scalarData;
    // scalarData.dt=dt;
    // scalarData.time=totalTime;
    // hsize_t dim[1];
    // dim[0]=sizeof(scalarData)/sizeof(exportValuesHDF5);
    // int rank=sizeof(dim)/sizeof(hsize_t);
    // CompType mtype(sizeof(exportValuesHDF5));
    // string test="dt";
    // mtype.insertMember(test,HOFFSET(exportValuesHDF5,dt),PredType::NATIVE_DOUBLE);
    // test ="time";
    // mtype.insertMember(test,HOFFSET(exportValuesHDF5,time),PredType::NATIVE_DOUBLE);
    // DataSpace space(rank,dim);
    // test=groupName+"/scalarData";
    // DataSet dataset = file.createDataSet(test,mtype,space);
    // dataset.write(&scalarData, mtype);

    file.close();
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
    error.printError();
    return -1;
   }

   return 0;  // successfully terminated
}

int multipleBeamSimulator::input_hdf5_restart(const std::string &fileName,const std::string &groupName)
{
   // Try block to detect exceptions raised by any of the calls inside it
  try
  {
    Exception::dontPrint();

    H5std_string FILE_NAME(fileName.c_str());
    H5File file(FILE_NAME,H5F_ACC_RDONLY);
    string dataName;
    string groupName1="/initial_setting";

    dataName=groupName+"/contact/u_stat_wall";
    HDF5input::read_3dDOUBLEARRAY(file,dataName,u_stat_wall);
    dataName=groupName+"/contact/wallContact";
    HDF5input::read_2dINTARRAY(file,dataName,wallContact);

    for(int ibeam=0;ibeam<numOfBeams;ibeam++){

      //read scalars

      //read arrays
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/ie";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].ie);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/R0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].R0);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/ql0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql0);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/l0";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].l0);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/Ma";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].Ma);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/Ia";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].Ia);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/x0";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x0);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/x_ref";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x_ref);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/t0";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t0);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/t_ref";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t_ref);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/externalForce";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalForce);
      dataName=groupName1+"/beam_"+to_string(ibeam)+"/externalTorque";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalTorque);

      dataName=groupName+"/beam_"+to_string(ibeam)+"/x";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].x);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/t";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].t);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/mask1";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].mask1);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/mask2";
      HDF5input::read_2dINTARRAY(file,dataName,beam[ibeam].mask2);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fl";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fl);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fg";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fg);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/force";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].f);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/Torque";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].T);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/fcon";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].fcon);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/u_local";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].u_local);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ql";
      HDF5input::read_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/ln";
      HDF5input::read_1dDOUBLEARRAY(file,dataName,beam[ibeam].ln);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/u";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].u);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/v";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].v);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/udot2";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].udot2);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/q";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].q);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/qv";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].qv);
      dataName=groupName+"/beam_"+to_string(ibeam)+"/qdot2";
      HDF5input::read_2dDOUBLEARRAY(file,dataName,beam[ibeam].qdot2);
      
    }
    // for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    //   groupName2=groupName+"/beam_" +to_string(ibeam);
    //   file.createGroup(groupName2.c_str());
    // }

    // //scalar value export
    // exportValuesHDF5 scalarData;
    // scalarData.dt=dt;
    // scalarData.time=totalTime;
    // hsize_t dim[1];
    // dim[0]=sizeof(scalarData)/sizeof(exportValuesHDF5);
    // int rank=sizeof(dim)/sizeof(hsize_t);
    // CompType mtype(sizeof(exportValuesHDF5));
    // string test="dt";
    // mtype.insertMember(test,HOFFSET(exportValuesHDF5,dt),PredType::NATIVE_DOUBLE);
    // test ="time";
    // mtype.insertMember(test,HOFFSET(exportValuesHDF5,time),PredType::NATIVE_DOUBLE);
    // DataSpace space(rank,dim);
    // test=groupName+"/scalarData";
    // DataSet dataset = file.createDataSet(test,mtype,space);
    // dataset.write(&scalarData, mtype);

    file.close();
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
    error.printError();
    return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
    error.printError();
    return -1;
   }

   return 0;  // successfully terminated
}

void multipleBeamSimulator::readBeamValuables(H5::H5File &file,const std::string &groupName,const int ibeam)
{
  exportBeamValuesHDF5 beamData;

  hsize_t dim[1];
  dim[0]=sizeof(beamData)/sizeof(exportBeamValuesHDF5);
  int rank=sizeof(dim)/sizeof(hsize_t);
  CompType mtype(sizeof(exportBeamValuesHDF5));

  string test="numOfNode";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,numOfNode),PredType::NATIVE_INT);
  test ="numOfElm";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,numOfElm),PredType::NATIVE_INT);
  test ="rho";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,rho),PredType::NATIVE_DOUBLE);
  test ="mu";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,mu),PredType::NATIVE_DOUBLE);
  test ="rad";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,rad),PredType::NATIVE_DOUBLE);
  test ="wireD";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,wireD),PredType::NATIVE_DOUBLE);
  test ="EA";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,EA),PredType::NATIVE_DOUBLE);
  test ="EI";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,EI),PredType::NATIVE_DOUBLE);
  test ="GJ";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,GJ),PredType::NATIVE_DOUBLE);
  test ="Young";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,Young),PredType::NATIVE_DOUBLE);
  test ="Poisson";
  mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,Poisson),PredType::NATIVE_DOUBLE);

  DataSpace space(rank,dim);
  test=groupName+"/beamValuables";
  DataSet dataset = file.createDataSet(test,mtype,space);
  dataset.read(&beamData, mtype);

  beam[ibeam].numOfNode=beamData.numOfNode;   //node number
  beam[ibeam].numOfElm=beamData.numOfElm;   //element number
  beam[ibeam].rho=beamData.rho; //density
  beam[ibeam].mu=beamData.mu;  //viscosity
  beam[ibeam].rad=beamData.rad; //radius of coil
  beam[ibeam].wireD=beamData.wireD; //wire diamter
  beam[ibeam].EA=beamData.EA;  //stiffness
  beam[ibeam].EI=beamData.EI;  //bending stiffness
  beam[ibeam].GJ=beamData.GJ;  //torsional stiffness
  beam[ibeam].Young=beamData.Young; //Young's modulus
  beam[ibeam].Poisson=beamData.Poisson; //Poisson's ratio

}
