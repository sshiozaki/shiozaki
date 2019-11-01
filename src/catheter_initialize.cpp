#include "multipleBeams.h"

using namespace std;

// void multipleBeamSimulator::initialize_catheter()
// {
//     set_parameters_catheter();

//     int i=numOfBeams;    //beam[numOfBeams].でカテーテル
//     beam[i].number=i;
//     beam[i].read_beam_geometry_catheter(tp,i);



//     //read_geometry_catheter(tp);

//     initializePhysicalValues_catheter();

//     inputParameters_catheter(tp);
// }

void multipleBeamSimulator::set_parameters_catheter()
{
    set_calculationCondition_catheter();
    set_outputCondition_catheter();
}

void multipleBeamSimulator::set_calculationCondition_catheter()
{
    string base_label, label;

    base_label = "/CalculationCondition_catheter";

    label = base_label + "/Restart_catheter";
    if ( !tp.getInspectedValue(label,Restart_catheter)){
        cout << label << " is not set" << endl;
        exit(0);
    }
    if(Restart_catheter=="yes"){
        label = base_label + "/RestartNumber_catheter";
        if ( !tp.getInspectedValue(label,restartNumber_catheter)){
        cout << label << " is not set" << endl;
        exit(0);
        }
        label = base_label + "/RestartFile_catheter";
        if ( !tp.getInspectedValue(label,restartFile_catheter)){
        cout << label << " is not set" << endl;
        exit(0);
        }
    }


    label = base_label + "/itermax_catheter";
    if(!tp.getInspectedValue(label, itermax_catheter)){
        cout << label << " is not set" << endl;
        exit(0);
    }

    label = base_label + "/DiameterofCatheter";
    if(!tp.getInspectedValue(label, D_catheter)){
        cout << label << " is not set" << endl;
        exit(0);
    }
 
}

void multipleBeamSimulator::set_outputCondition_catheter()
{
  string base_label,label;

  base_label = "/Output_catheter";
  label = base_label + "/outputDir_catheter";
  if ( !tp.getInspectedValue(label, outputDir_catheter)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  label = base_label + "/fileName_catheter";
  if ( !tp.getInspectedValue(label, fileName_catheter)){
    cout << label << " is not set" << endl;
    exit(0);
  }

}

