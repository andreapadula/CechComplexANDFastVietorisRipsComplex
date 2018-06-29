/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// include the local headers
#include                  <ttkCechComplex.h>
#include                  <ttkProgramBase.h>

int main(int argc, char **argv) {

  vtkProgram<ttkCechComplex> program;
  
  // TODO-1: 
  // specify local parameters to the TTK module with default values.
    bool option1 = false;
    bool option2 = false;
    bool option3 = false;
    double radius = 0.0;
  // end of TODO-1

  // TODO-2:
  // register these arguments to the command line parser
    program.parser_.setArgument("R", &radius,"Radius", true);
    program.parser_.setOption("O1", &option1,"Cell-1 to enable or disable");
    program.parser_.setOption("O2", &option2,"Cell-2 to enable or disable");
    program.parser_.setOption("O3", &option3,"Cell-3 to enable or disable");
  // end of TODO-2
  
  int ret = 0;
  ret = program.init(argc, argv);
 
  if(ret != 0)
    return ret;

  // TODO-3:
  // change here the arguments of the vtkWrapper that you want to update prior
  // to execution.
    program.ttkObject_->SetRadius(radius);
    program.ttkObject_->SetOption1(option1);
    program.ttkObject_->SetOption2(option2);
    program.ttkObject_->SetOption3(option3);
  // end of TODO-3
  
  // execute data processing
  ret = program.run();
  
  if(ret != 0)
    return ret;
 
  // save the output
  // optional TODO-4:
  // if you want a different kind of output, re-implement the function save().
  ret = program.save();
  /// end of optional TODO-4
  
  return ret;
}
