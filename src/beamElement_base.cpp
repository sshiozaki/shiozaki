#include "beamElement_base.h"
using namespace std;

// #################################################################
/**
 * @brief read beam info. from tp file
 * @param [in]     tp      TextParser file
 * @param [in]     number  beam number
 */
void beamElement_base::read_beam_geometry(TextParser &tp,const int number)
{
  string base="/Beam"+to_string(number);
  readGeometryFromFile(tp,base);
  set_initial_t(t,x);
  set_initial_t(t0,x0);
  set_initial_t(t_ref,x_ref);

  setBoundaryCondition(tp,base);
  set_physicalCondition(tp,base);
}

// #################################################################
/**
 * @brief read boundary condition setting
 * @param [in]     tp      TextParser file
 * @param [in]     base    beam number info.
 */
void beamElement_base::setBoundaryCondition(TextParser &tp,const std::string &base)
{
  string force,torque,mask_u,mask_q;
  string base_label,label;

  base_label = base+"/Boundary";

  label = base_label + "/force";
  if ( !tp.getInspectedValue(label,force)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/torque";
  if ( !tp.getInspectedValue(label,torque)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/mask_u";
  if ( !tp.getInspectedValue(label,mask_u)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/mask_q";
  if ( !tp.getInspectedValue(label,mask_q)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  fileIO::read_original_file_double(externalForce,numOfNode,3,force.c_str());
  fileIO::read_original_file_double(externalTorque,numOfNode,3,torque.c_str());
  fileIO::read_original_file_int(mask1,numOfNode,3,mask_u.c_str());
  fileIO::read_original_file_int(mask2,numOfNode,3,mask_q.c_str());
}

// #################################################################
/**
 * @brief read beam geometry
 * @param [in]     tp      TextParser file
 * @param [in]     base    beam number info.
 */
void beamElement_base::readGeometryFromFile(TextParser &tp,const std::string &base)
{
  string str,base_label,label,nodeFile,nodeReferenceFile,elementFile;

  base_label = base;
  label = base_label + "/nodeFile";
  if ( !tp.getInspectedValue(label,nodeFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/nodeReferenceFile";
  if ( !tp.getInspectedValue(label,nodeReferenceFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/elementFile";
  if ( !tp.getInspectedValue(label,elementFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  numOfNode = fileIO::CountNumbersOfTextLines(nodeFile);
  numOfElm = fileIO::CountNumbersOfTextLines(elementFile);

  tensor_initialize();

  fileIO::read_original_file_double(x,numOfNode,3,nodeFile);
  fileIO::read_original_file_double(x0,numOfNode,3,nodeFile);
  fileIO::read_original_file_double(x_ref,numOfNode,3,nodeReferenceFile);
  fileIO::read_original_file_int(ie,numOfElm,2,elementFile);

}

// #################################################################
/**
 * @brief initial t setting based on Crisfield, Non-linear finite element analysis of solid and structures, p.203 (1997).
 * @param [out]    t    unit basis assigned in node
 * @param [in]     x    node position vector
 */
void beamElement_base::set_initial_t(DOUBLEARRAY3 &t,const DOUBLEARRAY2 &x)
{
  double e01[3],e[3],v[3],dv,length,p1q1,p2q1,p3q1,tmp2,tmp3;
  e01[0] = 0e0; e01[1] = 0e0; e01[2] = 1e0;

  t[0][0][0]=1e0; t[0][0][1]=0e0; t[0][0][2]=0e0;
  t[0][1][0]=0e0; t[0][1][1]=1e0; t[0][1][2]=0e0;
  t[0][2][0]=0e0; t[0][2][1]=0e0; t[0][2][2]=1e0;

  for(int j=0;j<3;j++) e[j]=x[1][j]-x[0][j];
  length=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
  for(int j=0;j<3;j++) e[j]/=length;
  p1q1=e[0]*t[0][0][0]+e[1]*t[0][0][1]+e[2]*t[0][0][2];
  if(p1q1+1e0<1e-12){
    t[0][0][0]=-1e0; t[0][0][1]=0e0; t[0][0][2]=0e0;
    t[0][1][0]=0e0; t[0][1][1]=-1e0; t[0][1][2]=0e0;
    t[0][2][0]=0e0; t[0][2][1]=0e0; t[0][2][2]=1e0;
  }

  for(int ic=0;ic<numOfNode-1;ic++){

    for(int j=0;j<3;j++) e[j]=x[ic+1][j]-x[ic][j];
    length=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
    for(int j=0;j<3;j++) e[j]/=length;

    p1q1=e[0]*t[ic][0][0]+e[1]*t[ic][0][1]+e[2]*t[ic][0][2];
    p2q1=e[0]*t[ic][1][0]+e[1]*t[ic][1][1]+e[2]*t[ic][1][2];
    p3q1=e[0]*t[ic][2][0]+e[1]*t[ic][2][1]+e[2]*t[ic][2][2];

    for(int j=0;j<3;j++){
      t[ic][1][j]-=(t[ic][0][j]+e[j])*p2q1/(1e0+p1q1);
      t[ic][2][j]-=(t[ic][0][j]+e[j])*p3q1/(1e0+p1q1);
      t[ic][0][j]=e[j];
    }

    tmp2=t[ic][1][0]*t[ic][1][0]+t[ic][1][1]*t[ic][1][1]+t[ic][1][2]*t[ic][1][2];
    tmp3=t[ic][2][0]*t[ic][2][0]+t[ic][2][1]*t[ic][2][1]+t[ic][2][2]*t[ic][2][2];
    for(int j=0;j<3;j++){
      t[ic][1][j]/=tmp2;
      t[ic][2][j]/=tmp3;
    }

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) t[ic+1][i][j]=t[ic][i][j];
    }
  }
}
// #################################################################
/**
 * @brief allocation
 */
void beamElement_base::tensor_initialize()
{
  ie=Allocation::allocate2dINT(numOfElm,2);
  x=Allocation::allocate2dDOUBLE(numOfNode,3);
  x0=Allocation::allocate2dDOUBLE(numOfNode,3);
  x_ref=Allocation::allocate2dDOUBLE(numOfNode,3);
  t0=Allocation::allocate3dDOUBLE(numOfNode,3,3);
  t=Allocation::allocate3dDOUBLE(numOfNode,3,3);
  t_ref=Allocation::allocate3dDOUBLE(numOfNode,3,3);

  mask1=Allocation::allocate2dINT(numOfNode,3);
  mask2=Allocation::allocate2dINT(numOfNode,3);

  externalForce=Allocation::allocate2dDOUBLE(numOfNode,3);
  externalTorque=Allocation::allocate2dDOUBLE(numOfNode,3);

  //contact beam
  contact_beam2=Allocation::allocate2dINT(numOfElm,3);

  wire_cross=Allocation::allocate1dDOUBLE(numOfElm);
}

// #################################################################
/**
 * @brief set physical conditions
 * @param [in]     tp      TextParser file
 * @param [in]     base    beam number info.
 */
void beamElement_base::set_physicalCondition(TextParser &tp,const std::string &base)
{
  string base_label,label;

  base_label = base+"/PhysicalCondition";

  label = base_label + "/Young";
  if ( !tp.getInspectedValue(label,Young)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/G";
  if ( !tp.getInspectedValue(label,G)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  // label = base_label + "/Poisson";
  // if ( !tp.getInspectedValue(label,Poisson)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/EA";
  // if ( !tp.getInspectedValue(label,EA)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/EI";
  // if ( !tp.getInspectedValue(label,EI)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  // label = base_label + "/GJ";
  // if ( !tp.getInspectedValue(label,GJ)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }
  label = base_label + "/rho";
  if ( !tp.getInspectedValue(label,rho)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/mu";
  if ( !tp.getInspectedValue(label,mu)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/radius";
  if ( !tp.getInspectedValue(label,rad)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/wireDiameter";
  if ( !tp.getInspectedValue(label,wireD)){
    cout << label << " is not set" << endl;
    exit(0);
  }
}