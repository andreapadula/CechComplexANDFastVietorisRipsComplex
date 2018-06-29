/// \ingroup vtkWrappers
/// \class ttkFastVietorisRipsComplex
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the fastVietorisRipsComplex processing package.
///
/// VTK wrapping code for the @FastVietorisRipsComplex package.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::FastVietorisRipsComplex
#pragma once

// ttk code includes
#include                  <FastVietorisRipsComplex.h>
#include                  <ttkWrapper.h>

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
class VTKFILTERSCORE_EXPORT ttkFastVietorisRipsComplex 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
    
    static ttkFastVietorisRipsComplex* New();
    vtkTypeMacro(ttkFastVietorisRipsComplex, vtkDataSetAlgorithm)
    
    // default ttk setters
    vtkSetMacro(debugLevel_, int);
    
    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
      SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
        
    // TODO-4
    // set-getters macros to define from each variable you want to access from 
    // the outside (in particular from paraview) - to adapt.
    // Note that the XML file for the ParaView plug-in specification needs to be
    // edited accordingly.
      vtkSetMacro(Radius, double);
      vtkGetMacro(Radius, double);
      vtkSetMacro(Option1, bool);
      vtkGetMacro(Option1, bool);
      
      
      vtkSetMacro(Option2, bool);
      vtkGetMacro(Option2, bool);
      
      vtkSetMacro(Option3, bool);
      vtkGetMacro(Option3, bool);
    // end of TODO-4

    // TODO-2
    // Over-ride the input types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the input types, on a per input basis.
    // In this example, the first input type is forced to vtkUnstructuredGrid.
    // The second input type is forced to vtkImageData.
//     int FillInputPortInformation(int port, vtkInformation *info){
//       
//       switch(port){
//         case 0:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
//           break;
//         case 1:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
//           break;
//         default:
//           break;
//       }
//       
//       return 1;
//     }
    // end of TODO-2
    
    // TODO-3
    // Over-ride the output types.
    // By default, this filter has one input and one output, of the same type.
    // Here, you can re-define the output types, on a per output basis.
    // In this example, the first output type is forced to vtkUnstructuredGrid.
    // The second output type is forced to vtkImageData.
//     int FillOutputPortInformation(int port, vtkInformation *info){
//       
//       switch(port){
//         case 0:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid"); 
//           break;
//         case 1:
//           info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData"); 
//           break;
//         default:
//           break;
//       }
//       
//       return 1;
//     }
    // end of TODO-3
    
    
  protected:
   
    ttkFastVietorisRipsComplex(){
      
        // init
        Radius = -1;
        Option1=false;
        Option2=false;
        Option3=false;
        UseAllCores = false;
      
      // TODO-1
      // Specify the number of input and output ports.
      // By default, this filter has one input and one output.
      // In this example, we define 2 inputs and 2 outputs.
//       SetNumberOfInputPorts(2);
//       SetNumberOfOutputPorts(2);
      // end of TODO-1
    }
    
    ~ttkFastVietorisRipsComplex(){};
    
    TTK_SETUP();
    int FillOutputPortInformation(int port, vtkInformation* info);
    int FillInputPortInformation(int port, vtkInformation* info);
    void incrementalVR(const vector<vector<int> > &neighbors,int k);
      void addCofaces(const vector<vector<int> > &neighbors,int k, const vector<int> &index,const vector<int> &local);
      void intersection(const vector<int> &a,const vector<int> &b,vector<int> &result);

      
  private:
    
      double                   Radius;
      bool                    Option1;
      bool                    Option2;
      bool                    Option3;
      FastVietorisRipsComplex                 fastVietorisRipsComplex_;
      vector<vector<int>> triangles;
      vector<vector<int>> tetras;
      vector<vector<int>> edges;
    
};
