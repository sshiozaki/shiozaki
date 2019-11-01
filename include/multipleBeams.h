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
#ifndef _MULTIPLE_BEAMS_H_
#define _MULTIPLE_BEAMS_H_

/**
 * @file   multipleBeams.h
 * @brief  multiple beam treatment main Header
 * @author T. Otani
 */
#include "coRotationalBeam.h"
#include "wall.h"
#include <unordered_map>
#include "hdf5_input.h"
#include "hdf5_output.h"
#include "hdf5_dataStructure.h"
#include "vof_conv.h"

//for frictional contact
typedef struct{
  double u_static[3];
}ContactDisplacement;

//tangential contact history
class BeamContactHistory{
 public:
  std::unordered_map<int,ContactDisplacement> contactHistoryMAP;
};
// 
//contact data storage (temporal values)
typedef struct{
 public:
  int ic;
  bool stickSlip;  //0:slip 1:stick
  double fn[4][3],ft[4][3],u_static[3];
}BeamContactComponent;
//contact data storage (temporal values)
class BeamContactTMP{
 public:
  std::vector<BeamContactComponent> component;
};

class multipleBeamSimulator{
 public:
  multipleBeamSimulator(){
    dt=1e0;
    totalTime=0e0;
  }

  int OMPnumThreads;
  int numOfBeams;
  int itermax;


  double R;
  //int NW1, NW2;
  //int N;

  //int discritization_mesh;
  int release_node_number;
  //INTARRAY4 contact_beam;



  std::string fileName;
  std::string outputDir,restartDir;
  //restart setting
  int restartNumber;
  std::string restartFile,Restart;

  TextParser tp;

  coRotationalBeam *beam;
  Wall *wall;

  VOFC *VOFC;


  //catheter-----------------------------------------------------
  int itermax_catheter;
  double D_catheter;
  std::string fileName_catheter;
  std::string outputDir_catheter, restartDir_catheter;
  //restart setting
  int restartNumber_catheter;
  std::string restartFile_catheter,Restart_catheter;

  //int numOfNode_catheter;
  //DOUBLEARRAY2 x_catheter;
  //-------------------------------------------------------------



 private:
  double recordInterval;
  std::string aneurysmFile,balloonFile;
  double dt,totalTime;
  double recordTime;
  int itercount;

  double mu_static_aneurysm,mu_static_coil;
  double coefficientOfRestitution;
  double minimumSlidingVelocity;
  double stickFrictionViscosCoefficient;
  //Wall *wall;
  BeamContactHistory ***beamContactHistory;   //if number=3, components are [1,1],[1,2],[1,3],[2,2],[2,3],[3,3]. Little bit complicated...
  BeamContactTMP ***beamContactTMP;

 //multipleBeams.cpp
 public:
  //void vof_main();

  void initialize();
  void initialize_beam();
  void outputCurrentStatus();
  void preprocess();
  void preprocess_auto();
  void mainLoop_explicit();
  void mainLoop_explicitAPCScheme(const int count, const int countmin);
  void mainLoop_explicitAPCScheme_catheter(const int count, const int countmin);
  void mainLoop_quasiStaticScheme();

  void add_force_share_node(const int ibeam1);
  void wire_cross(const int ibeam1, double &wire_tmp, int &count);
  void calc_spring_force(const int ibeam1);
  void calc_spring_force_norm(const int ibeam1, const int ibeam2, const int ic);

  void exportPVTP_polyLine(const std::string &file, const int filenumber);
  void exportPVTP_polyLine_2(const std::string &file);

  void exportPVTP_polyLine_catheter(const std::string &file);

  //catheter------------------------------------------------------------------------
  void initialize_catheter();
  void set_parameters_catheter();
  void set_calculationCondition_catheter();
  void set_outputCondition_catheter();
  //void read_geometry_catheter(TextParser &tp);
  void readGeometryFromFile_catheter(TextParser &tp, std::string base);
  void initializePhysicalValues_catheter();
  void inputParameters_catheter(TextParser &tp);

  void calc_SDF_catheter();
  void calc_distance(const int ix, const int iy, const int iz);
  void make_catheter_surface();

  void exportVTIBinary(const std::string &file, const std::string &name);
  //--------------------------------------------------------------------------------

 private:

  void set_parameters();
  void set_outputCondition();
  void set_calculationCondition();
  void inputBasicParameters();
  double calcAneurysmVolume(Wall &aneurysm,Wall &balloon);

  //contactBox.cpp
  double origin[3];
  int numOfVoxel[3];
  double dx[3];
  INTARRAY3 box;
  INTARRAY4 numOfNodeInBox;
  INTARRAY5 nodeInBox;
  void initializeBox();
  void initializeBox_catheter();
  void register_box();
  void free_box();

  //wall_contact.cpp
  DOUBLEARRAY3 u_stat_wall,u_stat_wall_tmp;  //static friction (beam*node*dim)
  INTARRAY2 wallContact;  //contact judge for friction
  double contactDetection_SDF(double (&g)[3],int (&ix)[3],coRotationalBeam &beam,const DOUBLEARRAY2 &x,const int ic);
  double contactDetection_SDF_catheter(double (&g)[3],int (&ix)[3],coRotationalBeam &beam,const DOUBLEARRAY2 &x,const int ic);
  void calc_wallContactForce_SDF(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,const double kcon,const double vcon,const double (&normal)[3],const double length,const int ic);
  void calc_wallFrictionForce_SDF(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,const double E_astalisk,const double kcon,const double (&normal)[3],const double length,const int ibeam,const int ic);
  void calc_unitNormalVector(double (&normal)[3],const double (&g)[3],const int (&ix)[3]);
  void calc_unitNormalVector_catheter(double (&normal)[3],const double (&g)[3],const int (&ix)[3]);
  void corrector_wallStaticFriction();

  //wall_contact_catheter.cpp-------------------------------------------------
  DOUBLEARRAY2 u_stat_wall_catheter,u_stat_wall_tmp_catheter;  //static friction (beam*node*dim)
  INTARRAY1 wallContact_catheter;  //contact judge for friction

  //---------------------------------------------------------------------------

  //contact.cpp
  void closedSurfaceContact(const int ibeam);
  void closedSurfaceContact_APC(const int ibeam);
  void closedSurfaceContact_APC_catheter(const int ibeam);
  void closedSurfaceContact_APC_2(const int ibeam, const int release_node_number);
  void closedSurfaceContact_APC_catheter_2(const int ibeam, const int release_node_number);

  void selfContact(const int ibeam);
  void selfContact_APC(const int ibeam);

  void contact_beams(const int ibeam1,const int ibeam2);
  void contact_beams_APC(const int ibeam1,const int ibeam2);

  void corrector_beamStaticFriction();
  void ContactInBox(std::vector<BeamContactComponent> &component,const int ic1,const int ic2,coRotationalBeam &beam1,coRotationalBeam &beam2,std::unordered_map<int,ContactDisplacement> map,
  const double kcon,const double E_astalisk);
  void calc_wallContactForce(coRotationalBeam &beam,const DOUBLEARRAY2 &x,const DOUBLEARRAY2 &v,const int ic,const double Radius);

  void calc_beamContactForce(double (&fn)[4][3],const coRotationalBeam &beam1,const coRotationalBeam &beam2,const double s,const double t,const double kcon,const double vcon,const double distance,const double (&normal)[3],const double (&vn)[3]);
  void calc_beamContactForce_friction(double (&ft)[4][3],const coRotationalBeam &beam1,const coRotationalBeam &beam2,const double s,const double t,const double kcon,const double distance,
      const double (&normal)[3],const double (&vn)[3],double (&us)[3],const double E_astalisk,bool &stickSlip);
  double beam2beam_contact_search(const DOUBLEARRAY2 &x1,const DOUBLEARRAY2 &x2,
                coRotationalBeam &beam1,coRotationalBeam &beam2,const int ic1,const int ic2,double &s,double &t,double (&g)[3]);

  void search_contact_beam(const int ibeam1, const int ibeam2);
  void search_contact_beam_2(const int ibeam1, const int ic1);
  void wire_cross_force(const int ibeam1, const int ibeam2, const int ic1, const int ic2, const double exf_norm, const DOUBLEARRAY1 wire_cross);

  //output.cpp
 public:
  int input_hdf5_all(const std::string &fileName,const std::string &groupName);
  int input_hdf5_restart(const std::string &fileName,const std::string &groupName);
  int output_hdf5_valuables(const std::string &fileName,const std::string &groupName);
  int output_hdf5_constants(const std::string &fileName);

  //catheter------------------------------------------------------------------------------
  int input_hdf5_all_catheter(const std::string &fileName,const std::string &groupName);
  int input_hdf5_restart_catheter(const std::string &fileName,const std::string &groupName);
  //int output_hdf5_valuables_catheter(const std::string &fileName,const std::string &groupName);
  int output_hdf5_constants_catheter(const std::string &fileName);
  //--------------------------------------------------------------------------------------


 private:
  void readBeamValuables(H5::H5File &file,const std::string &groupName,const int ibeam);
  void writeBeamValuables(H5::H5File &file,const std::string &groupName,const int ibeam);
  void output_hdf5_beamContact(H5::H5File &file,const std::string &groupName);
};


// class catheterSimulator{
//     public:
//     //int OMPnumThreads;
//     int itermax_catheter;

//     double D_catheter;

//     //std::string fileName;
//     std::string outputDir_catheter;

//     TextParser tp;

//     public:
//     int numOfNode_catheter;
//     DOUBLEARRAY2 x_catheter;

//     public:
//     //double origin[3];
//     //double globalLength[3];
//     //int numOfVoxel[3];
//     //double dx[3];
//     //DOUBLEARRAY3 SDF;
//     //std::string SDFfile;


//     public:
//     void initialize_catheter();
//     void set_parameters_catheter();
//     void set_calculationCondition_catheter();
//     void set_outputCondition_catheter();
//     void read_geometry_catheter(TextParser &tp);
//     void readGeometryFromFile_catheter(TextParser &tp, std::string base);
//     void initializePhysicalValues_catheter();
//     void inputParameters_catheter(TextParser &tp);

//     void calc_SDF_catheter();
//     void calc_distance(const int ix, const int iy, const int iz);
//     void make_catheter_surface();

//     void exportVTIBinary(const std::string &file, const std::string &name);
// };



#endif //_beam_DEPLOYMENT_H_