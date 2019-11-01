/**
 * @file coRotationalBeam.cpp
 * @brief coRotationalBeam class
 * @author T. Otani
 */

#include "coRotationalBeam.h"
using namespace std;

// #################################################################
/**
 * @brief calc unit basis in all nodes
 * @param [out]     t      unit basis assigned in node
 * @param [in]      q1     rotation angle
 */
// void coRotationalBeam::calc_t(DOUBLEARRAY3 &t,const DOUBLEARRAY2 &q1)
// {
//   double tmp,R[3][3];

//   #pragma omp parallel for private(tmp,R)
//   for(int ic=0;ic<numOfNode;ic++){
//     calc_R(R,q1,ic);
//     for(int i=0;i<3;i++){

//       for(int j=0;j<3;j++){
//         t[ic][i][j] = 0e0;
//         for(int k=0;k<3;k++) t[ic][i][j] += R[j][k] * t0[ic][i][k];
//       }
//     tmp=sqrt(t[ic][i][0]*t[ic][i][0]+t[ic][i][1]*t[ic][i][1]+t[ic][i][2]*t[ic][i][2]);
//     if(abs(tmp-1e0)<1e-14) tmp=1e0;
//     for(int j=0;j<3;j++) t[ic][i][j]/=tmp;
//     }
//   }
// }

void coRotationalBeam::calc_t(DOUBLEARRAY3 &t,const DOUBLEARRAY3 &t_pre,const DOUBLEARRAY2 &q1)
{
  double tmp,R[3][3],vec[3];

  #pragma omp parallel for private(tmp,R,vec)
  for(int ic=0;ic<numOfNode;ic++){
    calc_R(R,q1,ic);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) vec[j]=t_pre[ic][i][j];
      for(int j=0;j<3;j++){
        t[ic][i][j] = 0e0;
        for(int k=0;k<3;k++) t[ic][i][j] += R[j][k] * vec[k];
      }
      tmp=sqrt(t[ic][i][0]*t[ic][i][0]+t[ic][i][1]*t[ic][i][1]+t[ic][i][2]*t[ic][i][2]);
      if(fabs(tmp-1e0)<1e-15) continue;
      for(int j=0;j<3;j++) t[ic][i][j]/=tmp;
    }
  }
}
// #################################################################
/**
 * @brief calc rotational matrix
 * @param [out]     R      Rotation matrix
 * @param [in]      q1     rotation angle
 * @param [in]      ic     node number
 */
void coRotationalBeam::calc_R(double (&R)[3][3],const DOUBLEARRAY2 &q1,const int &ic)
{
  double th1,th2,th3,qt[4],v[3],dv,qval,S[3][3],S2[3][3],I[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      I[i][j]=0e0;
    }
    I[i][i]=1e0;
  }

  qval=q1[ic][0]*q1[ic][0]+q1[ic][1]*q1[ic][1]+q1[ic][2]*q1[ic][2];

  if(qval<1e-20){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        R[i][j]=0e0;
      }
      R[i][i]=1e0;
    }
  }else{
    qval=sqrt(qval);
    qt[0]=cos(5e-1*qval);
    qt[1]=sin(5e-1*qval)*q1[ic][0]/qval;
    qt[2]=sin(5e-1*qval)*q1[ic][1]/qval;
    qt[3]=sin(5e-1*qval)*q1[ic][2]/qval;
    R[0][0]=qt[0]*qt[0]+qt[1]*qt[1]-5e-1;
    R[0][1]=qt[1]*qt[2]-qt[3]*qt[0];
    R[0][2]=qt[1]*qt[3]+qt[2]*qt[0];

    R[1][0]=qt[1]*qt[2]+qt[3]*qt[0];
    R[1][1]=qt[0]*qt[0]+qt[2]*qt[2]-5e-1;
    R[1][2]=qt[2]*qt[3]-qt[1]*qt[0];

    R[2][0]=qt[1]*qt[3]-qt[2]*qt[0];
    R[2][1]=qt[2]*qt[3]+qt[1]*qt[0];
    R[2][2]=qt[0]*qt[0]+qt[3]*qt[3]-5e-1;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        R[i][j]*=2e0;
      }
    }
  }
}
// #################################################################
/**
 * @brief predictor scheme of APC (Vetter et al., European Journal of Mechanics A/Solids, 2013)
 * @param [in]     dt      time increment
 */
void coRotationalBeam::predictor(const double dt)
{
  double beta=2.5e-1;
  double gamma=5e-1;

  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      u_p[ic][j] = u[ic][j]+dt*v[ic][j] + (double)mask1[ic][j] * (5e-1*dt*dt*(1e0-2e0*beta)*udot2[ic][j]);
      q_p[ic][j] = q[ic][j]+dt*qv[ic][j] + (double)mask2[ic][j] * (5e-1*dt*dt*(1e0-2e0*beta)*qdot2[ic][j]);
      v_p[ic][j] = v[ic][j] + (double)mask1[ic][j] * (dt*(1e0-gamma)*udot2[ic][j]);
      qv_p[ic][j] = qv[ic][j] + (double)mask2[ic][j] * (dt*(1e0-gamma)*qdot2[ic][j]);
      x_p[ic][j] = x0[ic][j]+u_p[ic][j];
    }
  }

  double tmp;
  double PI=acos(-1.0e0);
  #pragma omp parallel for private(tmp)
  for(int ic=0;ic<numOfNode;ic++){
    tmp=sqrt(q_p[ic][0]*q_p[ic][0]+q_p[ic][1]*q_p[ic][1]+q_p[ic][2]*q_p[ic][2]);
    if(tmp<2e0*PI) continue;
    for(int j=0;j<3;j++) q_p[ic][j] = (tmp-2e0*PI)*q_p[ic][j]/tmp;
  }

  calc_t(t_p,t0,q_p);

}
// #################################################################
/**
 * @brief corrector scheme of APC (Vetter et al., European Journal of Mechanics A/Solids, 2013)
 * @param [in]     dt      time increment
 */
void coRotationalBeam::corrector(const double dt)
{
  double beta=2.5e-1;
  double gamma=5e-1;

  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      u[ic][j] = u_p[ic][j]+dt*dt*beta*udot2_p[ic][j];
      q[ic][j] = q_p[ic][j]+dt*dt*beta*qdot2_p[ic][j];
      v[ic][j] = v_p[ic][j]+dt*gamma*udot2_p[ic][j];
      qv[ic][j] = qv_p[ic][j]+dt*gamma*qdot2_p[ic][j];
      x[ic][j] = x0[ic][j]+u[ic][j];
      udot2[ic][j] = udot2_p[ic][j];
      qdot2[ic][j] = qdot2_p[ic][j];
    }
  }

  double tmp;
  double PI=acos(-1.0e0);
  #pragma omp parallel for private(tmp)
  for(int ic=0;ic<numOfNode;ic++){
   tmp=sqrt(q[ic][0]*q[ic][0]+q[ic][1]*q[ic][1]+q[ic][2]*q[ic][2]);
   if(tmp<2e0*PI) continue;
   for(int j=0;j<3;j++) q[ic][j] = (tmp-2e0*PI)*q[ic][j]/tmp;
  }

  calc_t(t,t0,q);

}

// #################################################################
/**
 * @brief calc temporal accelaration (Vetter et al., European Journal of Mechanics A/Solids, 2013)
 * @param [in]     dt      time increment
 */
void coRotationalBeam::calc_udot2(const double dt)
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      udot2_p[ic][j] = (double)mask1[ic][j] * (-f[ic][j]+fcon[ic][j]+externalForce[ic][j]-mu*Ma[ic]*v_p[ic][j])/Ma[ic];
      qdot2_p[ic][j] = (double)mask2[ic][j] * (-T[ic][j]+externalTorque[ic][j]-mu*Ia[ic][j]*qv_p[ic][j])/Ia[ic][j];
    }
  }
}

// #################################################################
/**
 * @brief estimate error (Vetter et al., European Journal of Mechanics A/Solids, 2013)
 * @param [in]     uref     maximum error of displacement
 * @param [in]     qref     maximum error of rotational displacement
 * @param [in]     dt      time increment
 */
double coRotationalBeam::estimator(const double uref,const double qref,const double dt)
{
  double beta=2.5e-1;
  double eta,eta1,eta2,norm1,norm2,norm_max1=0e0,norm_max2=0e0;

  for(int ic=0;ic<numOfNode;ic++){
    norm1=0e0;
    norm2=0e0;
    for(int j=0;j<3;j++){
      norm1+=(udot2[ic][j]-udot2_p[ic][j])*(udot2[ic][j]-udot2_p[ic][j]);
      norm2+=(qdot2[ic][j]-qdot2_p[ic][j])*(qdot2[ic][j]-qdot2_p[ic][j]);
    }
    norm1=sqrt(norm1);
    norm2=sqrt(norm2);
    if(norm_max1<norm1) norm_max1=norm1;
    if(norm_max2<norm2) norm_max2=norm2;
  }

  eta1=(beta-1e0/6e0)*dt*dt/uref*norm_max1;
  eta2=(beta-1e0/6e0)*dt*dt/qref*norm_max2;

  if(eta1>eta2){
    eta=eta1;
  }else{
    eta=eta2;
  }

  return eta;
}

// #################################################################
/**
 * @brief calc quasi-static scheme
 * @param [in]     dt      time increment
 */
void coRotationalBeam::quasiStaticScheme(const double dt)
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      x[ic][j] += (double)mask1[ic][j] * (dt*(-f[ic][j]+externalForce[ic][j]))/mu;
      q[ic][j] += (double)mask2[ic][j] * (dt*(-T[ic][j]+externalTorque[ic][j]))/mu;
    }
  }
}

// #################################################################
/**
 * @brief calc semi-Implicit Euler scheme
 * @param [in]     dt      time increment
 */
void coRotationalBeam::semiImplicit_Euler_scheme(const double dt)
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      v[ic][j] = (double)mask1[ic][j] * (Ma[ic]*v[ic][j]+dt*(-f[ic][j]+fcon[ic][j]+externalForce[ic][j]))/(Ma[ic]+mu*Ma[ic]*dt);
      qv[ic][j] = (double)mask2[ic][j] * (Ia[ic][j]*qv[ic][j]+dt*(-T[ic][j]+externalTorque[ic][j]))/(Ia[ic][j]+mu*Ia[ic][j]*dt);
    }
    for(int j=0;j<3;j++){
      x[ic][j] += v[ic][j]*dt;
      q[ic][j] += qv[ic][j]*dt;
    }
  }
}

// #################################################################
/**
 * @brief calc element length at current coordinates
 * @param [in]     x      node position vector
 */
void coRotationalBeam::calc_current_length(const DOUBLEARRAY2 &x)
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    ln[ic] = sqrt(pow(x[ie[ic][1]][0]-x[ie[ic][0]][0],2e0)
          + pow(x[ie[ic][1]][1]-x[ie[ic][0]][1],2e0)
          + pow(x[ie[ic][1]][2]-x[ie[ic][0]][2],2e0));
  }
}
// #################################################################
/**
 * @brief calc inertia moment
 */
void coRotationalBeam::calc_current_inertia_moment()
{
  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0;i<3;i++) Ia[ic][i] = 0e0;
  }
}

// #################################################################
/**
 * @brief calc radius of all node 
 */
void coRotationalBeam::calc_r(DOUBLEARRAY1 &r, const DOUBLEARRAY2 &x)
{
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    r[i] = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
  }
}

// #################################################################
/**
 * @brief calc kinetic energy of all node 
 */
void coRotationalBeam::calc_kinetic_energy(DOUBLEARRAY1 &kene, const DOUBLEARRAY2 &v)
{
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    kene[i] = 5.0e-1*sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  }
}

// #################################################################
/**
 * @brief calc total kinetic energy  
 */
double coRotationalBeam::calc_total(DOUBLEARRAY1 &total1)
{
  double total=0e0;

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    total = total + total1[i];
  }  

  return total;
}

// #################################################################
/**
 * @brief calc max  
 */
double coRotationalBeam::calc_max(DOUBLEARRAY1 &max1, const double max2)
{
  double max=max2;

  for(int i=0;i<numOfNode;i++){
    if(max1[i] > max){
      max = max1[i];
    }
  }
  
  return max;
}


// #########################################################
/**
 * @brief export vtp file
 * @param [in]  file     vtp file name
 */
void coRotationalBeam::exportVTP_polyLine(const std::string &file)
{
  FILE *fp;
  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<PolyData>\n");
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",numOfNode,1);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");

  fprintf(fp,"<Lines>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\">\n");
  for(int i=0;i<numOfNode;i++) fprintf(fp,"%d ",i);
  fprintf(fp,"\n</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\">\n");
  fprintf(fp,"%d\n",numOfNode);
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</Lines>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",v[i][0],v[i][1],v[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"AngularVelocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",qv[i][0],qv[i][1],qv[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Angle\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",q[i][0],q[i][1],q[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",-T[i][0],-T[i][1],-T[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",-f[i][0],-f[i][1],-f[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"contactForce\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",fcon[i][0],fcon[i][1],fcon[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t[i][0][0],t[i][0][1],t[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t[i][1][0],t[i][1][1],t[i][1][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t3\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t[i][2][0],t[i][2][1],t[i][2][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"mask1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",mask1[i][0],mask1[i][1],mask1[i][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"mask2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",mask2[i][0],mask2[i][1],mask2[i][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>\n");  
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_longitudinal\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ic1=0; ic1<numOfElm; ic1++){
    fprintf(fp, "%e\n", energy_longitudinal[ic1]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_torsion\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ic1=0; ic1<numOfElm; ic1++){
    fprintf(fp, "%e\n", energy_torsion[ic1]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Energy_bending\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int ic1=0; ic1<numOfElm; ic1++){
    fprintf(fp, "%e\n", energy_bending[ic1]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</PolyData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

// #########################################################
/**
 * @brief export vtp file
 * @param [in]  file    vtp file name
 * @param [in]  x1      node position vector
 * @param [in]  t1      unit basis
 */
void coRotationalBeam::exportVTP_polyLine(const std::string &file,const DOUBLEARRAY2 &x1,const DOUBLEARRAY3 &t1)
{
  FILE *fp;
  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<PolyData>\n");
  fprintf(fp,"<Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n",numOfNode,1);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",x1[i][0],x1[i][1],x1[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");

  fprintf(fp,"<Lines>\n");
  fprintf(fp,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\">\n");
  for(int i=0;i<numOfNode;i++) fprintf(fp,"%d ",i);
  fprintf(fp,"\n</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\">\n");
  fprintf(fp,"%d\n",numOfNode);
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</Lines>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t1[i][0][0],t1[i][0][1],t1[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t1[i][1][0],t1[i][1][1],t1[i][1][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"t3\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",t1[i][2][0],t1[i][2][1],t1[i][2][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"mask1\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",mask1[i][0],mask1[i][1],mask1[i][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Int32\" Name=\"mask2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",mask2[i][0],mask2[i][1],mask2[i][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>\n");  
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</PolyData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}