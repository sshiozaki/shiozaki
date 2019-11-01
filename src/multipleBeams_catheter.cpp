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
void multipleBeamSimulator::mainLoop_explicitAPCScheme_catheter(const int count, const int countmin)
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

    file = outputDir + "/" + "Catheter_0"+".vtp";
    exportPVTP_polyLine_2(file);
    file = outputDir + "/" + "Catheter_state_0"+".vtp";
    exportPVTP_polyLine_2(file);
  }
  printf("%f\n", recordTime);

  string HDF5File=outputDir_catheter + "/" +fileName+".h5";
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


// ##
// #########################################################
/**
 * @brief export vtp file
 * @param [in]  file     vtp file name
 */
void multipleBeamSimulator::exportPVTP_polyLine_catheter(const std::string &file)
{
  FILE *fp;
  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<PolyData>\n");
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",beam[numOfBeams].numOfNode, beam[numOfBeams].numOfElm);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].x[ic1][0], beam[ibeam].x[ic1][1], beam[ibeam].x[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");

  fprintf(fp,"<Lines>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode-1; ic1++){
      fprintf(fp, "%d %d\n", ic1, ic1+1);
    }
  }
  fprintf(fp,"\n</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\">\n");
  for(int i=1; i<=(beam[numOfBeams].numOfNode-1)*(numOfBeams); i++) fprintf(fp, "%d\n", 2*i);
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</Lines>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
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
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].q[ic1][0], beam[ibeam].q[ic1][1], beam[ibeam].q[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", -beam[ibeam].T[ic1][0], -beam[ibeam].T[ic1][1], -beam[ibeam].T[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", -beam[ibeam].f[ic1][0], -beam[ibeam].f[ic1][1], -beam[ibeam].f[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].fcon[ic1][0], beam[ibeam].fcon[ic1][1], beam[ibeam].fcon[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall_normal\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].fn_wall[ic1][0], beam[ibeam].fn_wall[ic1][1], beam[ibeam].fn_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall_tangental\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].ft_wall[ic1][0], beam[ibeam].ft_wall[ic1][1], beam[ibeam].ft_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_wall\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].f_wall[ic1][0], beam[ibeam].f_wall[ic1][1], beam[ibeam].f_wall[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce_beams\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].f_contact[ic1][0], beam[ibeam].f_contact[ic1][1], beam[ibeam].f_contact[ic1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][0][0], beam[ibeam].t[ic1][0][1], beam[ibeam].t[ic1][0][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][1][0], beam[ibeam].t[ic1][1][1], beam[ibeam].t[ic1][1][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t3\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfNode; ic1++){
      fprintf(fp, "%e %e %e\n", beam[ibeam].t[ic1][2][0], beam[ibeam].t[ic1][2][1], beam[ibeam].t[ic1][2][2]);
    }
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>\n");  
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_longitudinal\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfElm; ic1++){
      fprintf(fp, "%e\n", beam[ibeam].energy_longitudinal[ic1]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_torsion\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
    for(int ic1=0; ic1<beam[ibeam].numOfElm; ic1++){
      fprintf(fp, "%e\n", beam[ibeam].energy_torsion[ic1]);
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_bending\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ibeam=numOfBeams; ibeam<numOfBeams+1; ibeam++){
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

