#include "FilterOperations.h"

void FilterOperations::binarize(int outside, int inside, 
               int threshLo, int threshHi,
               RegistrationEngine::SourceType::Pointer *src, 
               RegistrationEngine::SourceType::Pointer *dest)
{
  //
  BinaryFilterType::Pointer binFilt = BinaryFilterType::New();
                                  
  binFilt->SetOutsideValue(outside);
  binFilt->SetInsideValue(inside); 
  binFilt->SetLowerThreshold(threshLo);
  binFilt->SetUpperThreshold(threshHi);
  
  binFilt->SetInput((*src)->GetOutput());

  binFilt->Update();

  *dest = binFilt;
}
               
void FilterOperations::rescaleIntensity(int min, int max,
                      RegistrationEngine::SourceType::Pointer *src, 
                      RegistrationEngine::SourceType::Pointer *dest) 
{
  //
  RescalerType::Pointer rescaler = RescalerType::New();

  rescaler->SetInput((*src)->GetOutput());
  
  rescaler->SetOutputMinimum(min);
  rescaler->SetOutputMaximum(max);
  
  rescaler->Update();
  
  *dest = rescaler;
}
