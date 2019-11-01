//##################################################################################
//
// multiple beam simulator
//
// Copyright (c) 2018- Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
//                    All rights reserved.
//
//##################################################################################

/**
 * @file multipleBeams.cpp
 * @author T. Otani
 */
#include "multipleBeams.h"
#include "vof_conv.h"

using namespace std;

// #################################################################
/**
 * @brief explicit adaptive time stepping predictor-corrector (APC) scheme
 * @brief Vetter et al., European Journal of Mechanics A/Solids, 2013
 */
void multipleBeamSimulator::mainLoop_explicitAPCScheme(const int count, const int countmin)
{
  //double recordTime;
  double eta,eta_tmp,uref=1e10;
  string file;
  double qref=acos(-1.0e0)/8.0e0;

  double etaMin=1e-5;
  double etaMax=1e1*etaMin;
  double dtmax=1e-1;
  double etaBar=sqrt(etaMin*etaMax);

  double kene_total;
  double kene_max;
  double kene_max2;
  double kene_max_tmp=1e-5;
  double kene_max_threshold=1e5;

  int iter_record=0;

  if(count==countmin){
    recordTime=recordInterval;
    itercount=0;

    kene_total=0e0;
    kene_max=0e0;

    file = outputDir + "/" + "FDstent_0"+".vtp";
    exportPVTP_polyLine_2(file);
    file = outputDir + "/" + "FDstent_state_0"+".vtp";
    exportPVTP_polyLine_2(file);
  }
  printf("%f\n", recordTime);

  string HDF5File=outputDir + "/" +fileName+".h5";
  string groupName;

  printf("itercount=%d\n", itercount);

  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
  	if(uref>beam[ibeam].rad) uref=beam[ibeam].rad;
	}

  printf("kinetic_total= %e kinetic_max= %e\n\n", kene_total, kene_max);

  //---------main loop start--------------
  for(int iter=1;iter<=itermax;iter++){

  	while(10){

    	for(int ibeam=0;ibeam<numOfBeams;++ibeam){
  	  	beam[ibeam].predictor(dt);
  	  	beam[ibeam].calc_t(beam[ibeam].t_p,beam[ibeam].t0,beam[ibeam].q_p);
			}
    	for(int ibeam=0;ibeam<numOfBeams;++ibeam){
    		beam[ibeam].calc_elastic_force(beam[ibeam].x_p,beam[ibeam].t_p);

        #pragma omp parallel for
    		for(int i=0;i<beam[ibeam].numOfNode;i++){
    			for(int j=0;j<3;j++){
            beam[ibeam].fcon[i][j]=0e0;
            beam[ibeam].f_wall[i][j]=0e0;
            beam[ibeam].fn_wall[i][j]=0e0;
            beam[ibeam].ft_wall[i][j]=0e0;
            beam[ibeam].f_contact[i][j]=0e0;
          }
    		}
    	}

      free_box();
      register_box();

      for(int ibeam=0; ibeam<numOfBeams; ibeam++){
        closedSurfaceContact_APC(ibeam);
      }

    	for(int ibeam1=0;ibeam1<numOfBeams;++ibeam1){
	      	for(int ibeam2=ibeam1+1;ibeam2<numOfBeams;++ibeam2){
	        	 contact_beams_APC(ibeam1,ibeam2);
	      	}
    	}

      for(int ibeam=0; ibeam<numOfBeams; ibeam++){
        for(int i=0; i<beam[ibeam].numOfNode; i++){
          #pragma omp parallel for
          for(int j=0; j<3; j++){
            beam[ibeam].f_contact[i][j] = beam[ibeam].fcon[i][j] - beam[ibeam].f_wall[i][j];
          }
        }
      }

    	//calc_udot2
    	for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    		beam[ibeam].calc_udot2(dt);
		  }

		  eta=0e0;
    	for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    		eta_tmp=beam[ibeam].estimator(uref,qref,dt);
    		if(eta<eta_tmp) eta=eta_tmp;
		  }

    	if(eta>etaMax){
      	dt=pow(etaBar/eta,1e0/3e0)*dt;
    	}else if(eta<etaMin){
      	dt=pow(etaBar/eta,1e0/3e0)*dt;
      	if(dt>dtmax) dt=dtmax;
      	break;
    	}else{
    		break;
    	}
   }

    //corrector scheme for static Friction deflection
    corrector_wallStaticFriction();
    corrector_beamStaticFriction();

    kene_total = 0e0;
    kene_max = 0e0;
    for(int ibeam=0;ibeam<numOfBeams;++ibeam){
    	beam[ibeam].corrector(dt);
    	beam[ibeam].calc_t(beam[ibeam].t,beam[ibeam].t0,beam[ibeam].q);

      beam[ibeam].calc_kinetic_energy(beam[ibeam].kene,beam[ibeam].v);
      kene_total = kene_total + beam[ibeam].calc_total(beam[ibeam].kene);

      kene_max2 = kene_max; 
      kene_max = beam[ibeam].calc_max(beam[ibeam].kene,kene_max2);
		}

    if(iter < 5000){
      if(kene_max > kene_max_tmp){
        kene_max_tmp = kene_max;
        kene_max_threshold = kene_max*1e-2;
        if(kene_max_threshold < 1e-1){
          kene_max_threshold = 2e-1;
        }
      }
    }

    totalTime+=dt;
    itercount = itercount + 1;
    if(totalTime>recordTime){
      printf("%3d %9d %3d  iteration=%8d totalTime=%e dt= %e kinetic_thre= %e kinetic_max= %e\n", count, itercount, (int)((recordTime+1e-5)/recordInterval), iter, totalTime, dt, kene_max_threshold, kene_max);

      for(int ibeam=0;ibeam<numOfBeams;ibeam++){
        file = outputDir + "/" + fileName+"_"+to_string(beam[ibeam].number)+"_"+to_string((int)(recordTime/recordInterval))+".vtp";
        beam[ibeam].exportVTP_polyLine(file);
      }

      file = outputDir + "/" + "FDstent_"+to_string((int)(recordTime/recordInterval))+".vtp";
      exportPVTP_polyLine_2(file);
      recordTime+=recordInterval;

      if(kene_max < kene_max_threshold){
        printf("%3d %9d %3d  iteration=%8d totalTime=%e dt= %e kinetic_thre= %e kinetic_max= %e\n", count, itercount, (int)((recordTime+1e-5)/recordInterval), iter, totalTime, dt, kene_max_threshold, kene_max);
        file = outputDir + "/" + "FDstent_state_"+to_string(count)+".vtp";
        exportPVTP_polyLine_2(file);
        groupName=to_string(count);
        output_hdf5_valuables(HDF5File,groupName);
        break;
      }else if(iter == itermax){
        printf("%3d %9d %3d  iteration=%8d totalTime=%e dt= %e kinetic_thre= %e kinetic_max= %e\n", count, itercount, (int)((recordTime+1e-5)/recordInterval), iter, totalTime, dt, kene_max_threshold, kene_max);
        file = outputDir + "/" + "FDstent_state_"+to_string(count)+".vtp";
        exportPVTP_polyLine_2(file);
        groupName=to_string(count);
        output_hdf5_valuables(HDF5File,groupName);
        break;
      }
    }

    if(iter > 50000){
      if(kene_max < kene_max_threshold){
        printf("%3d %9d %3d  iteration=%8d totalTime=%e dt= %e kinetic_thre= %e kinetic_max= %e\n", count, itercount, (int)((recordTime+1e-5)/recordInterval), iter, totalTime, dt, kene_max_threshold, kene_max);
        file = outputDir + "/" + "FDstent_state_"+to_string(count)+".vtp";
        exportPVTP_polyLine_2(file);
        groupName=to_string(count);
        output_hdf5_valuables(HDF5File,groupName);
        break;
      }else if(iter == itermax){
        printf("%3d %9d %3d  iteration=%8d totalTime=%e dt= %e kinetic_thre= %e kinetic_max= %e\n", count, itercount, (int)((recordTime+1e-5)/recordInterval), iter, totalTime, dt, kene_max_threshold, kene_max);
        file = outputDir + "/" + "FDstent_state_"+to_string(count)+".vtp";
        exportPVTP_polyLine_2(file);
        groupName=to_string(count);
        output_hdf5_valuables(HDF5File,groupName);
        break;
      }  
    }
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void multipleBeamSimulator::calc_spring_force_norm(const int ibeam1, const int ibeam2, const int ic)
// {
//   double normal[3];
//   double normal1[3], normal2[3];
//   double distance;
  
//   double PI=acos(-1.0e0);
//   double kcon;
//   double E_astalisk;
//   double vcon;
//   double vn[3];
//   double fn;
//   double vtmp, vtmp1, vtmp2;
//   double vnormal, vnormal1, vnormal2;

//   double ks=1e-1;
//   double v_mu=1e-4;

//   for(int j=0; j<3; j++){
//     normal[j] = beam[ibeam2].x[ic][j] - beam[ibeam1].x[ic][j];
//   }
//   distance = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

//   for(int j=0; j<3; j++){
//     normal[j] = normal[j]/distance;
//     vn[j] = beam[ibeam2].v_p[j] - beam[ibeam1].v_p[j];
//   }

//   kcon=(1e0-beam[ibeam1].Poisson*beam[ibeam1].Poisson)/(beam[ibeam1].EA*PI*pow(beam[ibeam1].rad,2e0))
//        +(1e0-beam[ibeam2].Poisson*beam[ibeam2].Poisson)/(beam[ibeam2].EA*PI*pow(beam[ibeam2].rad,2e0));
  
//   E_astalisk=1e0/kcon;

//   kcon=2.5e-1*acos(-1.0e0)*E_astalisk;

//   vcon=-2.0e0*log(coefficientOfRestitution)*sqrt(beam[ibeam1].Ma[ic]*kcon/(PI*PI+log(coefficientOfRestitution)*log(coefficientOfRestitution)));

//   //fn = ks*fabs(distance - beam[ibeam1].ini_distance[ic]);
//   fn = ks*fabs(distance - 4e-2);
//   vtmp = vn[0]*normal[0]+vn[1]*normal[1]+vn[2]*normal[2];
//   vnormal = vcon*vtmp*v_mu;

//   if(distance > 4e-2){
//     for(int j=0; j<3; j++){
//       beam[ibeam1].spforce[ic][j] = fn*normal[j] + vnormal*normal[j];
//       beam[ibeam2].spforce[ic][j] = (-1e0)*fn*normal[j] + (-1e0)*vnormal*normal[j];
//     }          
//   }else{
//     for(int j=0; j<3; j++){
//       beam[ibeam1].spforce[ic][j] = 0e0;
//       beam[ibeam2].spforce[ic][j] = 0e0;
//     }
//   }

// }

// void multipleBeamSimulator::calc_spring_force(const int ibeam1)
// {
//   int ibeam2=ibeam1+NW1;

//   for(int ic=0; ic<beam[ibeam1].numOfNode; ic++){
//     if(ic%discritization_mesh == 0){
//       if(ic==0+discritization_mesh){
//         calc_spring_force_norm(ibeam1, ibeam2, ic);
//       }else if(ic==(beam[ibeam1].numOfNode-1-discritization_mesh)){
//         calc_spring_force_norm(ibeam1, ibeam2, ic);
//       }

//       ibeam2 += 1;
//       if(ibeam2 > numOfBeams-1){
//         ibeam2 = NW1;
//       }
//     }
//   }
// }



// void multipleBeamSimulator::add_force_share_node(const int ibeam1)
// {
//   int ibeam2=ibeam1+NW1;

//   for(int ic=0; ic<beam[ibeam1].numOfNode; ic++){
//     if(ic%discritization_mesh == 0){
//       for(int j=0; j<3; j++){
//         beam[ibeam1].f_add[ic][j] += beam[ibeam2].f[ic][j];
//         //beam[ibeam1].T_add[ic][j] += beam[ibeam2].T[ic][j];

//         beam[ibeam2].f_add[ic][j] += beam[ibeam1].f[ic][j];
//         //beam[ibeam2].T_add[ic][j] += beam[ibeam1].T[ic][j];

//         beam[ibeam2].fcon[ic][j] = beam[ibeam1].fcon[ic][j];
//       }
//       ibeam2 += 1;
//       if(ibeam2 > numOfBeams-1){
//         ibeam2 = NW1;
//       }
//     }
//   }
// }



// #########################################################
/**
 * @brief corrector scheme for wall static friction
 */
void multipleBeamSimulator::corrector_wallStaticFriction()
{
  #pragma omp parallel for
  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int ic=0;ic<beam[ibeam].numOfNode;ic++){
      if(wallContact[ibeam][ic]==0) continue;
      for(int j=0;j<3;j++) u_stat_wall[ibeam][ic][j]=u_stat_wall_tmp[ibeam][ic][j];
    }
  }
}

// #########################################################
/**
 * @brief corrector scheme for bram static friction
 */
void multipleBeamSimulator::corrector_beamStaticFriction()
{
  int icc;
  ContactDisplacement Tmp;

  #pragma omp parallel for private(Tmp,icc)
  for(int ibeam=0;ibeam<numOfBeams;ibeam++){
    for(int ibeam2=0;ibeam2<numOfBeams-ibeam;ibeam2++){                      // caution!!!!!!!!!!!!!!
      for(int ic1=0;ic1<beam[ibeam].numOfElm;ic1++){

        beamContactHistory[ibeam][ibeam2][ic1].contactHistoryMAP.clear();

        for(int ic2=0;ic2<beamContactTMP[ibeam][ibeam2][ic1].component.size();ic2++){

          if(beamContactTMP[ibeam][ibeam2][ic1].component[ic2].stickSlip==false) continue;

          icc=beamContactTMP[ibeam][ibeam2][ic1].component[ic2].ic;

          for(int j=0;j<3;j++) Tmp.u_static[j]=beamContactTMP[ibeam][ibeam2][ic1].component[ic2].u_static[j];
          beamContactHistory[ibeam][ibeam2][ic1].contactHistoryMAP[icc]=Tmp;
        }

      }
    }
  }
}

// #########################################################
/**
 * @brief wire cross search
 */
void multipleBeamSimulator::wire_cross(const int ibeam1, double &wire_tmp, int &count)
{
  if(ibeam1 != 0){
    beam[ibeam1].wire_cross[0] = beam[ibeam1-1].wire_cross[0] * (-1e0);
    wire_tmp = beam[ibeam1].wire_cross[0]; 
  }
  count = 0; 
  for(int ic1=1; ic1<beam[ibeam1].numOfElm; ic1++){
    if(beam[ibeam1].contact_beam2[ic1][0] != -1){
      count += 1;
      if(count%2 == 0){
        beam[ibeam1].wire_cross[ic1] = (-1e0)*wire_tmp;
        wire_tmp = beam[ibeam1].wire_cross[ic1];
      }else{
        beam[ibeam1].wire_cross[ic1] = wire_tmp;
        wire_tmp = beam[ibeam1].wire_cross[ic1];
      }
    }
  }
}


// #########################################################
/**
 * @brief initialize routine for multiple beam simulator
 */
void multipleBeamSimulator::initialize()
{
  set_parameters();
  beam = new coRotationalBeam[numOfBeams+1];  //beam[numOfBeams].でカテーテルを表す
  wall = new Wall[2];

  if(Restart=="yes"){
    for(int i=0;i<numOfBeams;i++){
      beam[i].number=i;
      beam[i].read_beam_geometry(tp,i);
      beam[i].initializePhysicalValues();
    }
    string groupName=to_string(restartNumber);
    input_hdf5_restart(restartFile,groupName);
  }else{
    for(int i=0;i<numOfBeams;i++){
      beam[i].number=i;
      beam[i].read_beam_geometry(tp,i);
      beam[i].initializePhysicalValues();
    }
  }

  //catheter-----------------------------------------------------------------
  printf("catheter initialize\n");
  set_parameters_catheter();
  
  int ic=numOfBeams;     //beam[numOfBeams].  catheter
  // beam[ic].number=ic;
  // beam[ic].read_beam_geometry(tp,ic);
  // beam[ic].initializePhysicalValues();

  if(Restart_catheter=="yes"){
    beam[ic].number=ic;
    beam[ic].read_beam_geometry(tp,ic);
    beam[ic].initializePhysicalValues();

    string groupName_catheter=to_string(restartNumber_catheter);
    input_hdf5_restart_catheter(restartFile_catheter,groupName_catheter);
  }else{
    beam[ic].number=ic;
    beam[ic].read_beam_geometry(tp,ic);
    beam[ic].initializePhysicalValues();    
  }

  //-------------------------------------------------------------------------

  wall[0].inputParameters(tp);               //wall[0] artery
  wall[1].inputParameters_catheter(tp);      //wall[1] catheter

  initializeBox();

  //contact_beam=Allocation::allocate4dINT(numOfBeams,numOfBeams,beam[0].numOfElm,2);

  //contact setting
  u_stat_wall=Allocation::allocate3dDOUBLE(numOfBeams,beam[0].numOfNode,3);
  u_stat_wall_tmp=Allocation::allocate3dDOUBLE(numOfBeams,beam[0].numOfNode,3);
  wallContact=Allocation::allocate2dINT(numOfBeams,beam[0].numOfNode);

  //contact setting for catheter----------------------------------------------
  ic=numOfBeams;
  u_stat_wall_catheter=Allocation::allocate2dDOUBLE(beam[ic].numOfNode,3);
  u_stat_wall_tmp_catheter=Allocation::allocate2dDOUBLE(beam[ic].numOfNode,3);
  wallContact_catheter=Allocation::allocate1dINT(beam[ic].numOfNode);
  //--------------------------------------------------------------------------


  //caution! Below matrixes are not square.
  beamContactHistory=new BeamContactHistory**[numOfBeams];
  for(int i=0;i<numOfBeams;i++){
    beamContactHistory[i]=new BeamContactHistory*[numOfBeams-i];
    for(int j=0;j<numOfBeams-i;j++){
      beamContactHistory[i][j]=new BeamContactHistory[beam[i].numOfElm];
    }
  }

  beamContactTMP=new BeamContactTMP**[numOfBeams];
  for(int i=0;i<numOfBeams;i++){
    beamContactTMP[i]=new BeamContactTMP*[numOfBeams-i];
    for(int j=0;j<numOfBeams-i;j++){
      beamContactTMP[i][j]=new BeamContactTMP[beam[i].numOfElm];
    }
  }




}

// #######################################a##########################
/**
 * @brief read parameters for co-rotational beam
 */
void multipleBeamSimulator::set_parameters()
{
  set_calculationCondition();
  set_outputCondition();
}

// #################################################################
/**
 * @brief set output conditions
 */
void multipleBeamSimulator::set_outputCondition()
{
  string base_label,label;

  base_label = "/Output";
  label = base_label + "/outputDir";
  if ( !tp.getInspectedValue(label, outputDir)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/fileName";
  if ( !tp.getInspectedValue(label, fileName)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/recordInterval";
  if ( !tp.getInspectedValue(label,recordInterval)){
    cout << label << " is not set" << endl;
    exit(0);
  }
}

// #################################################################
/**
 * @brief set calculation conditions
 */
void multipleBeamSimulator::set_calculationCondition()
{
  string base_label,label;

  base_label = "/CalculationCondition";
  label = base_label + "/OMPnumThreads";
  if ( !tp.getInspectedValue(label,OMPnumThreads)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  label = base_label + "/Restart";
  if ( !tp.getInspectedValue(label,Restart)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  if(Restart=="yes"){
    label = base_label + "/RestartNumber";
    if ( !tp.getInspectedValue(label,restartNumber)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/RestartFile";
    if ( !tp.getInspectedValue(label,restartFile)){
      cout << label << " is not set" << endl;
      exit(0);
    }
  }

  label = base_label + "/numberOfBeams";
  if ( !tp.getInspectedValue(label,numOfBeams)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/DiameterofStent";
  if ( !tp.getInspectedValue(label,R)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  // label = base_label + "/N";
  // if ( !tp.getInspectedValue(label,N)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/NW1";
  // if ( !tp.getInspectedValue(label,NW1)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/NW2";
  // if ( !tp.getInspectedValue(label,NW2)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/discritization_mesh";
  // if ( !tp.getInspectedValue(label,discritization_mesh)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  label = base_label + "/iterMax";
  if ( !tp.getInspectedValue(label,itermax)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/AneurysmFriction";
  if ( !tp.getInspectedValue(label,mu_static_aneurysm)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/CoilFriction";
  if ( !tp.getInspectedValue(label,mu_static_coil)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/CoefficientOfRestitution";
  if ( !tp.getInspectedValue(label,coefficientOfRestitution)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/MinimumSlidingVelocity";
  if ( !tp.getInspectedValue(label,minimumSlidingVelocity)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/StickFrictionViscosCoefficient";
  if ( !tp.getInspectedValue(label,stickFrictionViscosCoefficient)){
    cout << label << " is not set" << endl;
    exit(0);
  }
}

// #########################################################
/**
 * @brief export pvtp file
 * @param [in]  file     pvtp file name
 */
void multipleBeamSimulator::exportPVTP_polyLine(const std::string &file, const int filenumber)
{
  FILE *fp;
  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<PPolyData GhostLevel=\"%d\">\n", numOfBeams);
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",beam[0].numOfNode,1);
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"<PPoints>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"</PPoints>\n");

  fprintf(fp,"<PLines>\n");
  fprintf(fp,"<PDataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\">\n");
  fprintf(fp,"\n</PDataArray>\n");

  fprintf(fp,"<PDataArray type=\"Int32\" Name=\"offsets\" format=\"appended\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"</PLines>\n");
  fprintf(fp,"<PPointData>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"AngularVelocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"Angle\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"contactForce\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");
  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"t1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"t2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"<PDataArray type=\"Float64\" Name=\"t3\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"<PDataArray type=\"Int32\" Name=\"mask1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"<PDataArray type=\"Int32\" Name=\"mask2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  fprintf(fp,"</PDataArray>\n");

  fprintf(fp,"</PPointData>\n");

  fprintf(fp,"<PCellData>\n");  
  fprintf(fp,"</PCellData>\n");

  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    fprintf(fp, "<Piece Source=\"test_%d_%d.vtp\">\n", ibeam, filenumber);
  }
  fprintf(fp,"</PPolyData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

// #########################################################
/**
 * @brief export vtp file
 * @param [in]  file     vtp file name
 */
void multipleBeamSimulator::exportPVTP_polyLine_2(const std::string &file)
{
  FILE *fp;
  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<PolyData>\n");
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",beam[0].numOfNode*(numOfBeams), beam[0].numOfElm*(numOfBeams));
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].x[ic1][0], beam[ibeam].x[ic1][1], beam[ibeam].x[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");

  fprintf(fp,"<Lines>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode-1; ic1++){
      fprintf(fp, "%d %d\n", ibeam*beam[ibeam].numOfNode+ic1, ibeam*beam[ibeam].numOfNode+ic1+1);
    }
  }
  fprintf(fp,"\n</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\">\n");
  for(int i=1; i<=(beam[0].numOfElm)*(numOfBeams); i++) fprintf(fp, "%d\n", 2*i);
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</Lines>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].v[ic1][0], beam[ibeam].v[ic1][1], beam[ibeam].v[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"AngularVelocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].qv[ic1][0], beam[ibeam].qv[ic1][1], beam[ibeam].qv[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Angle\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].q[ic1][0], beam[ibeam].q[ic1][1], beam[ibeam].q[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", -beam[ibeam].T[ic1][0], -beam[ibeam].T[ic1][1], -beam[ibeam].T[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", -beam[ibeam].f[ic1][0], -beam[ibeam].f[ic1][1], -beam[ibeam].f[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].fcon[ic1][0], beam[ibeam].fcon[ic1][1], beam[ibeam].fcon[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall_normal\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].fn_wall[ic1][0], beam[ibeam].fn_wall[ic1][1], beam[ibeam].fn_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall_tangental\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].ft_wall[ic1][0], beam[ibeam].ft_wall[ic1][1], beam[ibeam].ft_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].f_wall[ic1][0], beam[ibeam].f_wall[ic1][1], beam[ibeam].f_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_beams\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].f_contact[ic1][0], beam[ibeam].f_contact[ic1][1], beam[ibeam].f_contact[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][0][0], beam[ibeam].t[ic1][0][1], beam[ibeam].t[ic1][0][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][1][0], beam[ibeam].t[ic1][1][1], beam[ibeam].t[ic1][1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t3\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][2][0], beam[ibeam].t[ic1][2][1], beam[ibeam].t[ic1][2][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>\n");  
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_longitudinal\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfElm; ic1++){
      fprintf(fp, "%e\n", beam[ibeam].energy_longitudinal[ic1]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_torsion\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfElm; ic1++){
      fprintf(fp, "%e\n", beam[ibeam].energy_torsion[ic1]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_bending\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=0; ibeam<numOfBeams; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfElm; ic1++){
      fprintf(fp, "%e\n", beam[ibeam].energy_bending[ic1]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</PolyData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

