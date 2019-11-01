
//##################################################################################
//
// Multiple co-Rotational beam elements simulator
//
// Copyright (c) 2017- Mechanical and Bioengineering Systems Lab.,
//                     Department of Mechanical Science and Bioengineering,
//                     Graduate School of Engineering Science,
//                     Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @author T. Otani
 */

#include "multipleBeams.h"
#include "vof_conv.h"
using namespace std;

int main(int argc,char *argv[]){

  multipleBeamSimulator MBS;
  VOFC VOFC;

  string output;

  if(argc!=2){
    cout << "Invalid input. Plase set tp file." << endl;
    return -1;
  }

  //read tp file
  std::string input_file = argv[1];
  int ierror;
  if ((ierror = MBS.tp.read(input_file)) != TP_NO_ERROR) {
    printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }

  //VOFC.vof_main();

  MBS.initialize();
  //MBS.initialize_catheter();


  omp_set_num_threads(MBS.OMPnumThreads);

  mkdir(MBS.outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
  string testFile=MBS.outputDir + "/" +MBS.fileName+".h5";
  string groupName="0";
  MBS.output_hdf5_constants(testFile);
  //MBS.output_hdf5_valuables(testFile,groupName);
  MBS.input_hdf5_all(testFile,groupName);

  //catheter------------------------------------------------------------------------------------------
  mkdir(MBS.outputDir_catheter.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
  string testFile_catheter=MBS.outputDir_catheter + "/" +MBS.fileName_catheter+".h5";
  string groupName_catheter="0";
  MBS.output_hdf5_constants_catheter(testFile_catheter);
  //MBS.output_hdf5_valuables(testFile,groupName);
  MBS.input_hdf5_all_catheter(testFile_catheter,groupName_catheter);
  //--------------------------------------------------------------------------------------------------

  cout << endl << "     main loop start     " <<endl << endl;

  for(int i=0;i<MBS.numOfBeams;i++){
    // output = MBS.outputDir + "/" +MBS.fileName+"_"+to_string(MBS.beam[i].number)+"_ref"+".vtp";
    // MBS.beam[i].exportVTP_polyLine(output,MBS.beam[i].x_ref,MBS.beam[i].t_ref);
    // cout << output << endl;
    output = MBS.outputDir + "/" +MBS.fileName+"_"+to_string(MBS.beam[i].number)+"_0"+".vtp";
    MBS.beam[i].exportVTP_polyLine(output);
    cout << output << endl;
  }

  // int countmin=0;
  // int countmax=211;
  // for(int count=countmin; count<countmax; count++){
  //   //VOFC.vof_main_2(count, countmax);
  //   //MBS.wall[0].inputParameters(MBS.tp);
  //   MBS.mainLoop_explicitAPCScheme_catheter(count, countmin);


  //   MBS.mainLoop_explicitAPCScheme(count, countmin);
  // }






  //MBS.mainLoop_explicitAPCScheme();

  return 0;
}
