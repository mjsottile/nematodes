/***********************************************************************
 * CommandIterationUpdate class: Called by optimizer to display the
 * state of the optimization while it occurs.
 ***********************************************************************/
#ifndef __COMMANDITERATIONUPDATE_H__
#define __COMMANDITERATIONUPDATE_H__

#include "itkCommand.h"
#include "itkRegularStepGradientDescentOptimizer.h"
 
class CommandIterationUpdate : public itk::Command {
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
  
protected:
  CommandIterationUpdate() {};
  
public:
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef const OptimizerType   *OptimizerPointer;
  
  void Execute(itk::Object *caller, const itk::EventObject & event) {
    Execute( (const itk::Object *)caller, event);
  }
  
  void Execute(const itk::Object * object, const itk::EventObject & event) {
    OptimizerPointer optimizer = 
    dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) ) {
      return;
    }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};

#endif /* __COMMANDITERATIONUPDATE_H__ */
