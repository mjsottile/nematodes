/*
 *  RigidWormRegister.h
 *
 *  Created by Matthew Sottile on 8/8/08.
 */

#ifndef __RIGIDWORMREGISTER_H__
#define __RIGIDWORMREGISTER_H__

#include "CommandIterationUpdate.h"

#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegistrationMethod.h"
#include "itkThresholdImageFilter.h"

#include "InputHandler.h"
#include "RegistrationEngine.h"

#include "sexp.h"
#include "sexp_ops.h"

using namespace std;

class RigidOptimizerParameters {
public:
  float scaleAngle;
  float scaleRCX, scaleRCY;
  float scaleTX, scaleTY;
  
  float maxStepLength, minStepLength;
  unsigned int iterationCount;

  void print() {
    cout << "=====================================================" << endl;
    cout << "RigidOptimizerParameters" << endl;
    cout << "------------------------" << endl;
    cout << "  ANGLE SCALE = " << scaleAngle << endl;
    cout << "  ROT SCALE= (" << scaleRCX << "," << scaleRCY << ")" << endl;
    cout << "  TRANS SCALE = (" << scaleTX << "," << scaleTY << ")" << endl;
    cout << "  MIN/MAX STEP = " << minStepLength << "/" << maxStepLength << 
      endl;
    cout << "  ITERS = " << iterationCount << endl;
    cout << "=====================================================" << endl;
  }
};

class RigidWormRegistration : public RegistrationEngine {
public:  
  typedef itk::CenteredRigid2DTransform <double>   TransformType;
  
  typedef itk::CenteredTransformInitializer <TransformType,ImageType,ImageType>
    TransformInitializerType;
  
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef OptimizerType::ScalesType                OptimizerScalesType;
  
  typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> 
    MetricType;
    
  typedef itk::LinearInterpolateImageFunction<ImageType, double>
    InterpolatorType;
    
  typedef itk::ImageRegistrationMethod<ImageType, ImageType>
    RegistrationType;

  typedef itk::ThresholdImageFilter<ImageType> ThresholdFilterType;

public:
  RigidWormRegistration();
  virtual ~RigidWormRegistration() { }
  
  void enableObserver();
  
  //
  // setup functions
  //
  void setOptimizationParams(RigidOptimizerParameters p);
  void setDefaultOptimizationParams();
  void setDefaultExtractors();  
  virtual void setupLRRegistration();
  virtual void setupCrossFrameRegistration();
  void setTransformParameters(TransformType::ParametersType p);

  void initializeByGeometry(bool flag);

  void extractParams(RigidOptimizerParameters *params, sexp_t *rp);

  void enableThreshFilter(PixelType loThresh, PixelType hiThresh);
  
  //
  // worker functions
  //
  virtual void doRegistration();
  
  //
  // post processing functions
  //
  TransformType::ParametersType getLastRegistrationTransform();
  
  void writeTransformedMoving(string lfname, string rfname,
    TransformType::ParametersType p, TransformType::ParametersType lrp);
  
  OptimizerType::MeasureType getOptimizerValue() { 
    return optimizer->GetValue();
  }
  
private:
  void getOneParam(sexp_t *sx, string pname, float *val, bool *hasParam);
    
  TransformType::ParametersType ComposeTransforms(
    TransformType::ParametersType a,
    TransformType::ParametersType b);
    
  void setupRegistration(SourceType::Pointer fixed,
    SourceType::Pointer moving);
    
  void writeTransformedFrame(string fname, TransformType::ParametersType p, 
    SourceType::Pointer src);
    
  void debug(string s) { cerr << "RWR DEBUG: " << s << endl; }

protected:
  MetricType::Pointer                 metric;
  OptimizerType::Pointer              optimizer;
  InterpolatorType::Pointer           interpolator;
  RegistrationType::Pointer           registration;
  TransformType::Pointer              transform;
  TransformInitializerType::Pointer   initializer;
  CommandIterationUpdate::Pointer     observer;
  ThresholdFilterType::Pointer        movingThreshFilter, fixedThreshFilter;
    
  SourceType::Pointer                 leftFixedSource, rightFixedSource;
  SourceType::Pointer                 leftMovingSource, rightMovingSource;

  RigidOptimizerParameters params;
  
  PixelType lowerThresh, highThresh;
  bool thresholdFilterEnabled;

  bool initByGeom;
  bool optParamsSet;
  bool sourcesSet;
  bool registrationSet;
};

#endif /* __RIGIDWORMREGISTER_H__ */
