#ifndef __FILTEROPERATIONS_H__
#define __FILTEROPERATIONS_H__

#include "RegistrationEngine.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

class FilterOperations {
 public:
  typedef RegistrationEngine::ImageType ImageType;
   
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> 
    RescalerType;
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> 
    BinaryFilterType;
  
 public:
  FilterOperations() { }
   
  void binarize(int outside, int inside, int threshLo, int threshHi,
                RegistrationEngine::SourceType::Pointer *src, 
                RegistrationEngine::SourceType::Pointer *dest);
                 
  void rescaleIntensity(int min, int max,
                        RegistrationEngine::SourceType::Pointer *src, 
                        RegistrationEngine::SourceType::Pointer *dest);
};

#endif /* __FILTEROPERATIONS_H__ */
