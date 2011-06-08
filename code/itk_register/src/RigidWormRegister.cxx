/*******************************************************************
 *  RigidWormRegister.cxx        
 *
 *  Created by Matthew Sottile on 8/8/08.
 *******************************************************************/

#include "RigidWormRegister.h"

/**
 * constructor
 */
RigidWormRegistration::RigidWormRegistration() 
{
  //
  // create a set of objects for the registration
  //
  metric = MetricType::New();
  optimizer = OptimizerType::New();
  interpolator = InterpolatorType::New();
  registration = RegistrationType::New();
  
  //
  // setup the helpers for the registration object
  //
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  
  //
  // create a new transformation and associate it with the registration
  // object
  //
  transform = TransformType::New();
  registration->SetTransform(transform);
  
  //
  // create a transform initializer.  the initializer is good to
  // try and get a good initial guess given the input images
  // geometry and mass.  far better than trying to set these manually.
  //
  initializer = TransformInitializerType::New();
  
  //
  // set flags to false that are used to make sure we don't do things like
  // the registration before specific operations have occurred, like setting
  // the sources up or setting parameters for the optimizer or transform
  //
  optParamsSet = false;
  sourcesSet = false;
  registrationSet = false;
  
  initByGeom = false;

  //
  // filtering
  //
  thresholdFilterEnabled = false;
  lowerThresh = 0;
  highThresh = 0;
}

/**
 * Set up an observer to watch the optimization and report the progress
 * to STDOUT.
 */
void
RigidWormRegistration::enableObserver() 
{
  observer = CommandIterationUpdate::New();
  optimizer->AddObserver(itk::IterationEvent(), observer);
}

/**
 * Setup threshold filtering
 */
void
RigidWormRegistration::enableThreshFilter(PixelType lo, PixelType hi)
{
  thresholdFilterEnabled = true;
  lowerThresh = lo;
  highThresh = hi;

  movingThreshFilter = ThresholdFilterType::New();
  movingThreshFilter->ThresholdBelow(lo);
  movingThreshFilter->ThresholdAbove(hi);
  movingThreshFilter->SetOutsideValue( lo + ((hi-lo)/2) );

  fixedThreshFilter = ThresholdFilterType::New();
  fixedThreshFilter->ThresholdBelow(lo);
  fixedThreshFilter->ThresholdAbove(hi);
  fixedThreshFilter->SetOutsideValue( lo + ((hi-lo)/2) );
}

/**
 * Populate the params class field with a set of defaults for the
 * transformation (all zeros = do nothing) and the scales, step bounds,
 * and max iteration count for the optimizer.  These are reasonable
 * defaults if you are not starting with any knowledge of a good starting
 * transform.
 */
void
RigidWormRegistration::setDefaultOptimizationParams() 
{
  params.scaleAngle = 1.0 / 100.0;
  params.scaleRCX = 1.0 / 300.0;
  params.scaleRCY = 1.0 / 300.0;
  params.scaleTX = 1.0 / 300.0;
  params.scaleTY = 1.0 / 300.0;
  
  params.minStepLength = 0.000001;
  params.maxStepLength = 1.0;
  params.iterationCount = 100;
  
  // set the flag saying we have a valid parameter set
  optParamsSet = true;
}

/**
 * Given a parameter set, set the class field to it.  Set the flag
 * indicating that we have a valid parameter set.
 */
void
RigidWormRegistration::setOptimizationParams(RigidOptimizerParameters p) 
{
  params = p;
  
  optParamsSet = true;
}

/**
 * Set up the default extractors.  This is assuming that we have two
 * images (moving and fixed) that have meaningful subimages on the
 * left and the right.  The actual regions corresponding to each
 * extractor are defined in the RegistrationEngine class.
 */
void
RigidWormRegistration::setDefaultExtractors() 
{
  leftFixedExtractor->SetInput(fixedImageReader->GetOutput());
  rightFixedExtractor->SetInput(fixedImageReader->GetOutput());
  leftMovingExtractor->SetInput(movingImageReader->GetOutput());
  rightMovingExtractor->SetInput(movingImageReader->GetOutput());
  
  // setup the sources.  these are used to know what the end of the
  // current pipeline is regardless of what actually is there.
  leftFixedSource = leftFixedExtractor;
  leftMovingSource = leftMovingExtractor;
  rightFixedSource = rightFixedExtractor;
  rightMovingSource= rightMovingExtractor;
  
  // set a flag saying we have valid sources.
  sourcesSet = true;
}

/**
 * Given two parameter sets a and b, create a new parameter set
 * representing T(a) composed with T(b).
 */
RigidWormRegistration::TransformType::ParametersType
RigidWormRegistration::ComposeTransforms(
  RigidWormRegistration::TransformType::ParametersType a,
  RigidWormRegistration::TransformType::ParametersType b)
{
  TransformType::ParametersType parameters;
  TransformType::Pointer t1 = TransformType::New();
  TransformType::Pointer t2 = TransformType::New();
  
  t1->SetParameters(a);
  t2->SetParameters(b);
  
  t1->Compose(t2,true);
  
  parameters = t1->GetParameters();
  
  return parameters;
}

/**
 * Write the transformed moving frame out using the given fixed->moving
 * transform parameter p and the left/right transform lr.  The left frame
 * is written with just transform p, and the right frame is written with
 * lr composed with p.
 */
void
RigidWormRegistration::writeTransformedMoving(string lfname, string rfname,
		    RigidWormRegistration::TransformType::ParametersType p, 
		    RigidWormRegistration::TransformType::ParametersType lr)
{
  // Cast the rightMovingExtractor to a SourceType pointer
  SourceType::Pointer src = (SourceType::Pointer)rightMovingExtractor;

  // Force an update to read the frame.
  src->Update();
  
  //
  // Get the spacing from src.  This is necessary to know how to
  // scale the rightX and rightY offsets.  Cannot assume scale is 1!
  // This most often occurs with TIFF images - not with PNG or JPG.
  //
  ImageType::SpacingType spacing = 
    src->GetOutput()->GetSpacing();

  /*** NOTE: Left/Right transform disabled.
       TODO: make this an option in the input file
   ***/  
  //TransformType::ParametersType lr_mod = lr;
  //lr_mod[3] -= rightXOffset*spacing[0];
  //TransformType::ParametersType composed = ComposeTransforms(lr_mod,p);

  // Since L/R disabled, just make composed the input transform p
  TransformType::ParametersType composed = p;

  // Apply scaled offsets
  composed[1] += rightXOffset*spacing[0];
  composed[2] += rightYOffset*spacing[1];

  // Write left frame
  writeTransformedFrame(lfname,p,(SourceType::Pointer)leftMovingExtractor);
  
  // Write right frame
  writeTransformedFrame(rfname,composed,
			(SourceType::Pointer)rightMovingExtractor);
}

/**
 * Given a general fixed and moving source (end of the processing pipeline
 * for each), set the image inputs for the registration method, update
 * the pipelines, and set the initializer up.  This is called before the
 * do registration routine, and is required to allow us to reuse the
 * registration logic for both the Left/Right registration task and the
 * frame->frame task.
 */
void
RigidWormRegistration::setupRegistration(SourceType::Pointer fixed,
                                         SourceType::Pointer moving) 
{
  if (thresholdFilterEnabled) {
    fixedThreshFilter->SetInput(fixed->GetOutput());
    movingThreshFilter->SetInput(moving->GetOutput());
    registration->SetFixedImage(fixedThreshFilter->GetOutput());
    registration->SetMovingImage(movingThreshFilter->GetOutput());
  } else {
    registration->SetFixedImage(fixed->GetOutput());
    registration->SetMovingImage(moving->GetOutput());
  }

  fixed->Update();
  moving->Update();
  
  registration->SetFixedImageRegion(fixed->GetOutput()->GetBufferedRegion());
  initializer->SetTransform(transform);
  initializer->SetFixedImage(fixed->GetOutput());
  initializer->SetMovingImage(moving->GetOutput());
  if (initByGeom)
    initializer->GeometryOn();
  else
    initializer->MomentsOn();
  initializer->InitializeTransform();
}

void
RigidWormRegistration::initializeByGeometry(bool flag) {
  initByGeom = flag;
}

/**
 * Sets up the registration for Left/Right registration.
 */
void 
RigidWormRegistration::setupLRRegistration() 
{
  if (!sourcesSet) {
    debug("SOURCES NOT SET YET!");
    exit(EXIT_FAILURE);
  }
  setupRegistration(leftFixedSource, rightFixedSource);

  registrationSet = true;
}

/**
 * Sets up the registration for cross frame registration.
 */
void 
RigidWormRegistration::setupCrossFrameRegistration() 
{
  if (!sourcesSet) {
    debug("SOURCES NOT SET YET!");
    exit(EXIT_FAILURE);
  }
  setupRegistration(leftFixedSource, leftMovingSource);

  registrationSet = true;
}

/**
 * perform the registration
 */
void
RigidWormRegistration::doRegistration() 
{
  if (!optParamsSet) {
    debug("SET PARAMETERS FIRST!");
    exit(EXIT_FAILURE);
  }
  
  if (!sourcesSet) {
    debug("SET SOURCES FIRST!");
    exit(EXIT_FAILURE);
  }
  
  if (!registrationSet) {
	  debug("SET REGISTRATION UP FIRST!");
	  exit(EXIT_FAILURE);
  }
  
  registration->SetInitialTransformParameters(transform->GetParameters());  
  OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
    
  optimizerScales[0] = params.scaleAngle;
  optimizerScales[1] = params.scaleRCX;
  optimizerScales[2] = params.scaleRCY;
  optimizerScales[3] = params.scaleTX;
  optimizerScales[4] = params.scaleTY;
  
  optimizer->SetScales(optimizerScales);
  optimizer->SetMaximumStepLength(params.maxStepLength);
  optimizer->SetMinimumStepLength(params.minStepLength);
  optimizer->SetNumberOfIterations(params.iterationCount);

  fixedImageReader->Update();
  movingImageReader->Update();

  try {
    registration->StartRegistration();
  }  catch (itk::ExceptionObject &err) {
    debug("itk::ExceptionObject caught!");
    cerr << err << endl;
    exit(EXIT_FAILURE);
  }
}

RigidWormRegistration::TransformType::ParametersType
RigidWormRegistration::getLastRegistrationTransform()
{
  return registration->GetLastTransformParameters();
}

void
RigidWormRegistration::setTransformParameters(
  RigidWormRegistration::TransformType::ParametersType p)
{
  transform->SetParameters(p);
}

void 
RigidWormRegistration::getOneParam(
  sexp_t *sx, 
  string pname,
  float *val, 
  bool *hasParam) 
{
  sexp_t *param;
  char pstr[1024];
  string s;
  
  param = find_sexp(pname.c_str(), sx);
  
  if (param == NULL) {
    cout << "PARAM \"" << pname << "\" NOT FOUND." << endl;
    *hasParam = false;
  }
  
  *hasParam = true;
  
  if (param->next->ty == SEXP_VALUE) {
    strncpy(pstr,param->next->val,param->next->val_used);
    pstr[param->next->val_used] = '\0';
    *val = (float)strtof(pstr,NULL);
  } else {
    cerr << "INVALID PARAMETER VALUE: " << pname << endl;
    exit(EXIT_FAILURE);
  }
}
  
void 
RigidWormRegistration::extractParams(RigidOptimizerParameters *params, 
  sexp_t *rp)
{
  float val;
  bool hasParam;

  getOneParam(rp, "angle", &val, &hasParam);
  if (hasParam)
    params->scaleAngle = val;
  
  getOneParam(rp, "rotCenterX", &val, &hasParam);
  if (hasParam)
    params->scaleRCX = val;

  getOneParam(rp, "rotCenterY", &val, &hasParam);
  if (hasParam) 
    params->scaleRCY = val;

  getOneParam(rp, "translateX", &val, &hasParam);
  if (hasParam)
    params->scaleTX = val;

  getOneParam(rp, "translateY", &val, &hasParam);
  if (hasParam)
    params->scaleTY = val;

  getOneParam(rp, "maxStepLength", &val, &hasParam);
  if (hasParam)
    params->maxStepLength = val;
  
  getOneParam(rp, "minStepLength", &val, &hasParam);
  if (hasParam) 
    params->minStepLength = val;

  getOneParam(rp, "iterationCount", &val, &hasParam);
  if (hasParam) 
    params->iterationCount = (int)val;
}

/**
 * Given a filename, a set of transform parameters, and an image source,
 * write out a file to the filename with the image after the transform is
 * applied.
 */
void
RigidWormRegistration::writeTransformedFrame(string fname, 
					     RigidWormRegistration::TransformType::ParametersType p, 
					     SourceType::Pointer src) 
{
  //
  // make a new transform object
  //
  TransformType::Pointer xform = TransformType::New();
  
  //
  // force an update on the source
  //
  src->Update();
  
  //
  // inform the transform object of the new params
  //
  xform->SetParameters(p);
  
  //
  // make a resampler
  //
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  
  //
  // give the resampler the transform and the image source
  //
  resampler->SetTransform(xform);
  resampler->SetInput(src->GetOutput());
    
  //
  // get a pointer to the image itself from the source
  //
  ImageType::Pointer img = src->GetOutput();
  
  //
  // setup the resampler
  //
  resampler->SetSize(img->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(img->GetOrigin());
  resampler->SetOutputSpacing(img->GetSpacing());
  //  resampler->SetDefaultPixelValue(255);
  
  //
  // write the image output from the resampler
  //
  writeImage(fname,(SourceType::Pointer)resampler);
}
