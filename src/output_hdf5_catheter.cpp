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
 * @file output_hdf5.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <H5Cpp.h>

using namespace std;
using namespace H5;

// #########################################################
/**
 * @brief output constant values (beam constants)
 */
// int multipleBeamSimulator::output_hdf5_valuables(const std::string &fileName,const std::string &groupName)
// {
//    // Try block to detect exceptions raised by any of the calls inside it
//   try
//   {
//     Exception::dontPrint();

//     H5std_string FILE_NAME(fileName.c_str());
//     H5File file(FILE_NAME,H5F_ACC_RDWR);
//     Group group(file.createGroup(groupName.c_str()));
//     // file.openGroup(groupName.c_str());
//     string dataName;
//     string groupName2;

//     groupName2=groupName+"/contact";
//     file.createGroup(groupName2.c_str());
//     dataName=groupName2+"/u_stat_wall";
//     HDF5output::write_3dDOUBLEARRAY(file,dataName,u_stat_wall,numOfBeams,beam[0].numOfNode,3);
//     dataName=groupName2+"/wallContact";
//     HDF5output::write_2dINTARRAY(file,dataName,wallContact,numOfBeams,beam[0].numOfNode);

//     //beam contact values
//     output_hdf5_beamContact(file,groupName2);

//     for(int ibeam=0;ibeam<numOfBeams;ibeam++){
//       groupName2=groupName+"/beam_" +to_string(ibeam);
//       file.createGroup(groupName2.c_str());
//     }

//     for(int ibeam=0;ibeam<numOfBeams;ibeam++){
//       groupName2=groupName+"/beam_" +to_string(ibeam);
//       file.openGroup(groupName2.c_str());

//       dataName=groupName2+"/x";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].x,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/t";
//       HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].t,beam[ibeam].numOfNode,3,3);
//       dataName=groupName2+"/mask1";
//       HDF5output::write_2dINTARRAY(file,dataName,beam[ibeam].mask1,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/mask2";
//       HDF5output::write_2dINTARRAY(file,dataName,beam[ibeam].mask2,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/fl";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].fl,beam[ibeam].numOfNode,7);
//       dataName=groupName2+"/fg";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].fg,beam[ibeam].numOfNode,12);
//       dataName=groupName2+"/force";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].f,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/Torque";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].T,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/fcon";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].fcon,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/u_local";
//       HDF5output::write_1dDOUBLEARRAY(file,dataName,beam[ibeam].u_local,beam[ibeam].numOfElm);
//       dataName=groupName2+"/ql";
//       HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql,beam[ibeam].numOfElm,2,3);
//       dataName=groupName2+"/ln";
//       HDF5output::write_1dDOUBLEARRAY(file,dataName,beam[ibeam].ln,beam[ibeam].numOfElm);
//       dataName=groupName2+"/u";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].u,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/v";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].v,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/udot2";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].udot2,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/q";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].q,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/qv";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].qv,beam[ibeam].numOfNode,3);
//       dataName=groupName2+"/qdot2";
//       HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].qdot2,beam[ibeam].numOfNode,3);

//     }

//     //scalar value export
//     exportValuesHDF5 scalarData;
//     scalarData.dt=dt;
//     scalarData.time=totalTime;

//     hsize_t dim[1];
//     dim[0]=sizeof(scalarData)/sizeof(exportValuesHDF5);
//     int rank=sizeof(dim)/sizeof(hsize_t);
//     CompType mtype(sizeof(exportValuesHDF5));
//     string test="dt";
//     mtype.insertMember(test,HOFFSET(exportValuesHDF5,dt),PredType::NATIVE_DOUBLE);
//     test ="time";
//     mtype.insertMember(test,HOFFSET(exportValuesHDF5,time),PredType::NATIVE_DOUBLE);
//     DataSpace space(rank,dim);
//     test=groupName+"/scalarData";
//     DataSet dataset = file.createDataSet(test,mtype,space);
//     dataset.write(&scalarData, mtype);

//     file.close();
//    }  // end of try block

//    // catch failure caused by the H5File operations
//    catch( FileIException error )
//    {
//     error.printError();
//     return -1;
//    }

//    // catch failure caused by the DataSet operations
//    catch( DataSetIException error )
//    {
//     error.printError();
//     return -1;
//    }

//    // catch failure caused by the DataSpace operations
//    catch( DataSpaceIException error )
//    {
//     error.printError();
//     return -1;
//    }

//    // catch failure caused by the DataSpace operations
//    catch( DataTypeIException error )
//    {
//     error.printError();
//     return -1;
//    }

//    return 0;  // successfully terminated
// }

// // #########################################################
// /**
//  * @brief output beam contact values
//  */
// void multipleBeamSimulator::output_hdf5_beamContact(H5::H5File &file,const std::string &groupName)
// {
//   //scalar value export
//   exportBeamContactHDF5 *scalarData;
//   int number,rank;
//   hsize_t dim[1];

//   CompType mtype(sizeof(exportBeamContactHDF5));
//   string test="ic1";
//   mtype.insertMember(test,HOFFSET(exportBeamContactHDF5,ic1),PredType::NATIVE_INT);
//   test ="ic2";
//   mtype.insertMember(test,HOFFSET(exportBeamContactHDF5,ic2),PredType::NATIVE_INT);
//   test ="ux";
//   mtype.insertMember(test,HOFFSET(exportBeamContactHDF5,ux),PredType::NATIVE_DOUBLE);
//   test ="uy";
//   mtype.insertMember(test,HOFFSET(exportBeamContactHDF5,uy),PredType::NATIVE_DOUBLE);
//   test ="uz";
//   mtype.insertMember(test,HOFFSET(exportBeamContactHDF5,uz),PredType::NATIVE_DOUBLE);

//   for(int ibeam=0;ibeam<numOfBeams;ibeam++){
//     for(int ibeam2=0;ibeam2<numOfBeams-ibeam;ibeam2++){
//       number=0e0;
//       for(int i=0;i<beam[ibeam].numOfElm;i++) number+=beamContactHistory[ibeam][ibeam2][i].contactHistoryMAP.size();
//       scalarData = new exportBeamContactHDF5[number];

//       number=0e0;
//       for(int i=0;i<beam[ibeam].numOfElm;i++){
//         scalarData[number].ic1=i;
//         for(auto itr=beamContactHistory[ibeam][ibeam2][i].contactHistoryMAP.begin();itr!=beamContactHistory[ibeam][ibeam2][i].contactHistoryMAP.end();++itr){
//           scalarData[number].ic2=itr->first;
//           scalarData[number].ux=itr->second.u_static[0];
//           scalarData[number].uy=itr->second.u_static[1];
//           scalarData[number].uz=itr->second.u_static[2];
//         }
//       }
//       dim[0]=sizeof(scalarData)/sizeof(exportBeamContactHDF5);
//       rank=sizeof(dim)/sizeof(hsize_t);
//       DataSpace space(rank,dim);
//       test=groupName+"/"+to_string(ibeam)+"_"+to_string(ibeam2);
//       DataSet dataset(file.createDataSet(test,mtype,space));
//       dataset.write(scalarData,mtype);
//       number++;
//       delete scalarData;
//     }
//   }

// }

// #########################################################
/**
 * @brief output constant values
 */
int multipleBeamSimulator::output_hdf5_constants_catheter(const std::string &fileName)
{
    string dataName;
    string groupName="/initial_setting";

   // Try block to detect exceptions raised by any of the calls inside it
  try
  {
    Exception::dontPrint();

    H5std_string FILE_NAME(fileName.c_str());
    H5File file(FILE_NAME,H5F_ACC_TRUNC);
    Group group(file.createGroup(groupName.c_str()));
    // file.openGroup(groupName.c_str());

    for(int ibeam=numOfBeams;ibeam<numOfBeams+1;ibeam++){
      groupName="/initial_setting/beam_" +to_string(ibeam);
      file.createGroup(groupName.c_str());
    }

    for(int ibeam=numOfBeams;ibeam<numOfBeams+1;ibeam++){
      groupName="/initial_setting/beam_" +to_string(ibeam);
      file.openGroup(groupName.c_str());

      writeBeamValuables(file,groupName,ibeam);

      dataName=groupName+"/ie";
      HDF5output::write_2dINTARRAY(file,dataName,beam[ibeam].ie,beam[ibeam].numOfElm,2);
      dataName=groupName+"/R0";
      HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].R0,beam[ibeam].numOfElm,3,3);
      dataName=groupName+"/ql0";
      HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].ql0,beam[ibeam].numOfElm,2,3);
      dataName=groupName+"/l0";
      HDF5output::write_1dDOUBLEARRAY(file,dataName,beam[ibeam].l0,beam[ibeam].numOfElm);
      dataName=groupName+"/Ma";
      HDF5output::write_1dDOUBLEARRAY(file,dataName,beam[ibeam].Ma,beam[ibeam].numOfNode);
      dataName=groupName+"/Ia";
      HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].Ia,beam[ibeam].numOfNode,3);
      dataName=groupName+"/x0";
      HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].x0,beam[ibeam].numOfNode,3);
      dataName=groupName+"/x_ref";
      HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].x_ref,beam[ibeam].numOfNode,3);
      dataName=groupName+"/t0";
      HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].t0,beam[ibeam].numOfNode,3,3);
      dataName=groupName+"/t_ref";
      HDF5output::write_3dDOUBLEARRAY(file,dataName,beam[ibeam].t_ref,beam[ibeam].numOfNode,3,3);
      dataName=groupName+"/externalForce";
      HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalForce,beam[ibeam].numOfNode,3);
      dataName=groupName+"/externalTorque";
      HDF5output::write_2dDOUBLEARRAY(file,dataName,beam[ibeam].externalTorque,beam[ibeam].numOfNode,3);
      dataName=groupName+"/mask1";
      HDF5output::write_2dINTARRAY(file,dataName,beam[ibeam].mask1,beam[ibeam].numOfNode,3);
      dataName=groupName+"/mask2";
      HDF5output::write_2dINTARRAY(file,dataName,beam[ibeam].mask2,beam[ibeam].numOfNode,3);
    }

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

// #########################################################
/**
 * @brief output constant values (beam constants)
 */
// void multipleBeamSimulator::writeBeamValuables(H5::H5File &file,const std::string &groupName,const int ibeam)
// {
//   exportBeamValuesHDF5 beamData;

//   beamData.numOfNode=beam[ibeam].numOfNode;   //node number
//   beamData.numOfElm=beam[ibeam].numOfElm;   //element number
//   beamData.rho=beam[ibeam].rho; //density
//   beamData.mu=beam[ibeam].mu;  //viscosity
//   beamData.rad=beam[ibeam].rad; //radius of coil
//   beamData.wireD=beam[ibeam].wireD; //wire diamter
//   beamData.EA=beam[ibeam].EA;  //stiffness
//   beamData.EI=beam[ibeam].EI;  //bending stiffness
//   beamData.GJ=beam[ibeam].GJ;  //torsional stiffness
//   beamData.Young=beam[ibeam].Young; //Young's modulus
//   beamData.Poisson=beam[ibeam].Poisson; //Poisson's ratio

//   hsize_t dim[1];
//   dim[0]=sizeof(beamData)/sizeof(exportBeamValuesHDF5);
//   int rank=sizeof(dim)/sizeof(hsize_t);
//   CompType mtype(sizeof(exportBeamValuesHDF5));

//   string test="numOfNode";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,numOfNode),PredType::NATIVE_INT);
//   test ="numOfElm";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,numOfElm),PredType::NATIVE_INT);
//   test ="rho";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,rho),PredType::NATIVE_DOUBLE);
//   test ="mu";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,mu),PredType::NATIVE_DOUBLE);
//   test ="rad";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,rad),PredType::NATIVE_DOUBLE);
//   test ="wireD";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,wireD),PredType::NATIVE_DOUBLE);
//   test ="EA";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,EA),PredType::NATIVE_DOUBLE);
//   test ="EI";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,EI),PredType::NATIVE_DOUBLE);
//   test ="GJ";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,GJ),PredType::NATIVE_DOUBLE);
//   test ="Young";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,Young),PredType::NATIVE_DOUBLE);
//   test ="Poisson";
//   mtype.insertMember(test,HOFFSET(exportBeamValuesHDF5,Poisson),PredType::NATIVE_DOUBLE);

//   DataSpace space(rank,dim);
//   test=groupName+"/beamValuables";
//   DataSet dataset = file.createDataSet(test,mtype,space);
//   dataset.write(&beamData, mtype);
// }
