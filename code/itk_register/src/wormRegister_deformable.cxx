/***************************************************************/
/* wormRegister.cxx : An ITK-based tool for registering worm   */
/* movies from Shawn Lockery.                                  */
/*                                                             */
/* Matt Sottile / matt@cs.uoregon.edu                          */
/***************************************************************/

#include "itkTimeProbesCollectorBase.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkCommand.h"
#include "itkConstantPadImageFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkLBFGSOptimizer.h"

#include "CommandIterationUpdate.h"

using namespace std;

void usage(char *av0) {
  cout << "usage: " << av0 << " -t IMTY -b BASE -i DIR -o DIR -r FR# -s FR# -e FR#" << endl;
  cout << "   IMTY : Image type (aka: image extension)" << endl;
  cout << "   BASE : Filename base (name minus frame number)" << endl;
  cout << endl;
  exit(EXIT_FAILURE);
}

string buildFilename(string dir, string base, int num, string ext) {
  char fileno[10];
  string filename;
  sprintf(fileno,"%04d",num);
  filename = dir + "/" + base + string(fileno) + "." + ext;
  return filename;
}

/***********************************************************************
 * main()
 ***********************************************************************/
int main(int argc, char **argv) {
  //------------------------------------------------------------
  // parameter setup
  //------------------------------------------------------------
  
  // basic parameters
  int referenceframe = -1;
  int startframe = -1;
  int endframe = -1;
  string inDir = "";
  string outDir = "";
  string fnameBase = "";
  string fnameExt = "";
  
  // command line processing
  int opt = 0;
  opt = getopt(argc, argv, "t:b:r:s:e:i:o:?h");
  while (opt != -1) {	
    switch(opt) {
      case 'r':
        referenceframe = atoi(optarg);
        break;
      case 's':
        startframe = atoi(optarg);
        break;
      case 'e':
        endframe = atoi(optarg);
        break;
      case 't':
        fnameExt = string(optarg);
        break;
      case 'i':
        inDir = string(optarg);
        break;
      case 'o':
        outDir = string(optarg);
        break;
      case 'b':
        fnameBase = string(optarg);
        break;
      case 'h':
      case '?':
        usage(argv[0]);
        break;
    }
    opt = getopt(argc, argv, "t:b:r:s:e:i:o:?h");
  }
  
  if (referenceframe == -1 || startframe == -1 || endframe == -1) {
    usage(argv[0]);
  }
  
  // setup dimension and pixel type
  const    unsigned int    Dimension = 2;
  typedef  unsigned short    PixelType;
  
  // Define image types
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  typedef itk::Image< PixelType, Dimension >  ImageType;

  // setup file readers
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  
  FixedImageReaderType::Pointer  fixedImageReader  = 
    FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = 
    MovingImageReaderType::New();
  
  // set files
  string filename = buildFilename(inDir,fnameBase,referenceframe,fnameExt);
  fixedImageReader->SetFileName(  filename );
  fixedImageReader->Update();
  
  //=================================================================
  //=================================================================
  // Region of interest definition
  //=================================================================
  //=================================================================

  // set up a region of interest filter to extract left and right subimages
  typedef itk::RegionOfInterestImageFilter<ImageType,
                                           ImageType>    ExtractorType;

  // get the size of the frames from the fixed image reader
  ImageType::SizeType imgSize = 
    fixedImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();

  // Compute sizes of both halves of the image.  index 0 = horiz,
  // index 1 = vert

  //
  // sizes of just middle third (horizontal slice)
  //

  /*
  ImageType::IndexType leftStart;
  leftStart[0] = 0; leftStart[1] = imgSize[1]/3;
  ImageType::SizeType leftSize;
  leftSize[0] = imgSize[0]/2; leftSize[1] = imgSize[1]/3;
  
  ImageType::IndexType rightStart;
  rightStart[0] = leftStart[0]+1; rightStart[1] = imgSize[1]/3;
  ImageType::SizeType rightSize;
  rightSize[0] = imgSize[0]/2; rightSize[1] = imgSize[1]/3;
  */

  //
  // sizes for whole half images
  //
  ImageType::IndexType leftWholeStart;
  leftWholeStart[0] = 0; leftWholeStart[1] = 0;
  ImageType::SizeType leftWholeSize;
  leftWholeSize[0] = imgSize[0]/2; leftWholeSize[1] = imgSize[1];
  
  ImageType::IndexType rightWholeStart;
  rightWholeStart[0] = leftWholeStart[0]+1; rightWholeStart[1] = 0;
  ImageType::SizeType rightWholeSize;
  rightWholeSize[0] = imgSize[0]/2; rightWholeSize[1] = imgSize[1];

  // setup regions for windowed-third sub images
  /*
  ImageType::RegionType leftRegion, rightRegion;
  leftRegion.SetSize(leftSize);
  leftRegion.SetIndex(leftStart);
  rightRegion.SetSize(rightSize);
  rightRegion.SetIndex(rightStart);
  */

  // setup regions for windowed-third sub images
  ImageType::RegionType leftWholeRegion, rightWholeRegion;
  leftWholeRegion.SetSize(leftWholeSize);
  leftWholeRegion.SetIndex(leftWholeStart);
  rightWholeRegion.SetSize(rightWholeSize);
  rightWholeRegion.SetIndex(rightWholeStart);


  //
  // setup extractors for horizontal slices
  //
  // left
  /*
  ExtractorType::Pointer leftFixedExtractor = ExtractorType::New();
  leftFixedExtractor->SetRegionOfInterest(leftRegion);

  ExtractorType::Pointer leftMovingExtractor = ExtractorType::New();
  leftMovingExtractor->SetRegionOfInterest(leftRegion);

  // right
  ExtractorType::Pointer rightFixedExtractor = ExtractorType::New();
  rightFixedExtractor->SetRegionOfInterest(rightRegion);
  
  ExtractorType::Pointer rightMovingExtractor = ExtractorType::New();
  rightMovingExtractor->SetRegionOfInterest(rightRegion);
  */
  //
  // setup extractors for whole halves
  //
  // left
  ExtractorType::Pointer leftWholeFixedExtractor = ExtractorType::New();
  leftWholeFixedExtractor->SetRegionOfInterest(leftWholeRegion);

  ExtractorType::Pointer leftWholeMovingExtractor = ExtractorType::New();
  leftWholeMovingExtractor->SetRegionOfInterest(leftWholeRegion);

  // right
  ExtractorType::Pointer rightWholeFixedExtractor = ExtractorType::New();
  rightWholeFixedExtractor->SetRegionOfInterest(rightWholeRegion);
  
  ExtractorType::Pointer rightWholeMovingExtractor = ExtractorType::New();
  rightWholeMovingExtractor->SetRegionOfInterest(rightWholeRegion);
  

  //=================================================================
  //=================================================================
  // Registration pipeline setup
  //=================================================================
  //=================================================================

  // define transform, optimizer, metric, interpolator, and registration
  
  //
  // CenteredRigid2DTransform
  //
  typedef itk::CenteredRigid2DTransform< double > LRTransformType;
  
  //
  // RegularStepGradientDescentOptimizer
  //
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef OptimizerType::ScalesType    OptimizerScalesType;

  //
  // MeanSquaresImageToImageMetric
  //
  typedef itk::MeanSquaresImageToImageMetric< FixedImageType, 
                                              MovingImageType >
    MetricType;
    
  //
  // LinearInterpolateImageFunction
  //                                                   
  typedef itk::LinearInterpolateImageFunction< MovingImageType, 
                                               double >
    InterpolatorType;
                                                       
  typedef itk::ImageRegistrationMethod< FixedImageType, 
                                        MovingImageType >
                                                       RegistrationType;
  
  // create instances
  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   lrRegistration  = RegistrationType::New();
  
  // set up registration
  lrRegistration->SetMetric(        metric        );
  lrRegistration->SetOptimizer(     optimizer     );
  lrRegistration->SetInterpolator(  interpolator  );
  
  // create transform
  LRTransformType::Pointer  transform = LRTransformType::New();
  lrRegistration->SetTransform( transform );

  // variables for transform results.  initialize them here to zeros

  // first set of variables: frame to frame transformation
  double frameRotationCenterX = 0.0;
  double frameRotationCenterY = 0.0;
  double frameTranslationX = 0.0;
  double frameTranslationY = 0.0;
  double frameAngle = 0.0;

  // second set of variables: left to right transformation.  this
  // is STATIC and determined once based on the fixed image frame.
  double lrRotationCenterX = 0.0;
  double lrRotationCenterY = 0.0;
  double lrTranslationX = 0.0;
  double lrTranslationY = 0.0;
  double lrAngle = 0.0;
    
  // transform initializer
  typedef itk::CenteredTransformInitializer< LRTransformType, 
                                             FixedImageType, 
                                             MovingImageType >  
    TransformInitializerType;
  
  TransformInitializerType::Pointer initializer = 
    TransformInitializerType::New();
  
  // add observer to watch optimization
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // preprocessing filter
  //  typedef itk::GradientMagnitudeImageFilter<FixedImageType,FixedImageType>
  //  FixedGradFilterType;
  //typedef itk::GradientMagnitudeImageFilter<MovingImageType,MovingImageType>
  //  MovingGradFilterType;
  
  //FixedGradFilterType::Pointer fixedPreFilter = FixedGradFilterType::New();
  //MovingGradFilterType::Pointer movingPreFilter = MovingGradFilterType::New();

  //fixedPreFilter->SetInput(fixedImageReader->GetOutput());

  /*
  leftFixedExtractor->SetInput(fixedImageReader->GetOutput());
  leftFixedExtractor->Update();
  rightFixedExtractor->SetInput(fixedImageReader->GetOutput());
  rightFixedExtractor->Update();
  */

  // set up extractors for whole left/right frames
  leftWholeFixedExtractor->SetInput(fixedImageReader->GetOutput());
  leftWholeFixedExtractor->Update();
  rightWholeFixedExtractor->SetInput(fixedImageReader->GetOutput());
  rightWholeFixedExtractor->Update();

  //=================================================================
  //=================================================================
  // Perform preliminary registration of left to right side
  //=================================================================
  //=================================================================
  
/*
  typedef itk::ConstantPadImageFilter<FixedImageType, FixedImageType> 
    PadType;
  PadType::Pointer leftPad = PadType::New();
  PadType::Pointer rightPad = PadType::New();
*/
    
  // first, create a RescaleIntensityImageFilter to scale the intensities
  // of both sides.

  typedef itk::RescaleIntensityImageFilter <FixedImageType, 
                                            FixedImageType> RescalerType;

  RescalerType::Pointer leftIntensityRescaler = RescalerType::New();
  RescalerType::Pointer rightIntensityRescaler = RescalerType::New();

  leftIntensityRescaler->SetInput(leftWholeFixedExtractor->GetOutput());
  leftIntensityRescaler->SetOutputMinimum(0);
  leftIntensityRescaler->SetOutputMaximum(255);

  rightIntensityRescaler->SetInput(rightWholeFixedExtractor->GetOutput());
  rightIntensityRescaler->SetOutputMinimum(0);
  rightIntensityRescaler->SetOutputMaximum(255);

/*
  leftPad->SetPadLowerBound( &leftWholeSize[0] );
  leftPad->SetPadUpperBound( &leftWholeSize[0] );
  leftPad->SetConstant( 0 );
  leftPad->SetInput( leftWholeFixedExtractor->GetOutput() );
//  leftPad->SetInput( leftIntensityRescaler->GetOutput() );

  rightPad->SetPadLowerBound( &leftWholeSize[0] );
  rightPad->SetPadUpperBound( &leftWholeSize[0] );
  rightPad->SetConstant( 0 );
  rightPad->SetInput( rightWholeFixedExtractor->GetOutput() );
//  rightPad->SetInput( rightIntensityRescaler->GetOutput() );
*/

//  registration->SetFixedImage(leftPad->GetOutput());
//  registration->SetMovingImage(rightPad->GetOutput());
  lrRegistration->SetFixedImage(leftIntensityRescaler->GetOutput());
  lrRegistration->SetMovingImage(rightIntensityRescaler->GetOutput());
  leftIntensityRescaler->Update();
  rightIntensityRescaler->Update();  
//  leftPad->Update();
//  rightPad->Update();

  lrRegistration->SetFixedImageRegion(
     	leftIntensityRescaler->GetOutput()->GetBufferedRegion() );
	
  initializer->SetTransform( transform );
  initializer->SetFixedImage(  leftIntensityRescaler->GetOutput() );
  initializer->SetMovingImage( rightIntensityRescaler->GetOutput() );
  
  /*
  // binarize the image to perform registration based only on the zero
  // values outside the image
  typedef itk::BinaryThresholdImageFilter< 
     FixedImageType, FixedImageType > BinaryFilterType;
     
  BinaryFilterType::Pointer leftBinaryFilter = BinaryFilterType::New();
  BinaryFilterType::Pointer rightBinaryFilter = BinaryFilterType::New();
  
  leftBinaryFilter->SetOutsideValue( 255 ); 
  leftBinaryFilter->SetInsideValue( 1 ); 
  leftBinaryFilter->SetLowerThreshold( 0 ); 
  leftBinaryFilter->SetUpperThreshold( 10 );
  
  rightBinaryFilter->SetOutsideValue( 255 ); 
  rightBinaryFilter->SetInsideValue( 1 ); 
  rightBinaryFilter->SetLowerThreshold( 0 ); 
  rightBinaryFilter->SetUpperThreshold( 10 );
  
  leftBinaryFilter->SetInput(leftIntensityRescaler->GetOutput());
  rightBinaryFilter->SetInput(rightIntensityRescaler->GetOutput());
  registration->SetFixedImage(leftBinaryFilter->GetOutput());
  registration->SetMovingImage(rightBinaryFilter->GetOutput());
  leftBinaryFilter->Update();
  rightBinaryFilter->Update();  

  registration->SetFixedImageRegion(
	   leftBinaryFilter->GetOutput()->GetBufferedRegion() );
	
  initializer->SetTransform( transform );
  initializer->SetFixedImage(  leftBinaryFilter->GetOutput() );
  initializer->SetMovingImage( rightBinaryFilter->GetOutput() );
  */      
  initializer->MomentsOn();
  
  initializer->InitializeTransform();
  transform->SetAngle( 0.0 );
  
  lrRegistration->SetInitialTransformParameters( transform->GetParameters() );

  // setup optimizer
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  const double lrTranslationScale = 1.0 / 300.0;

  optimizerScales[0] =  1.0 / 10.0;
  optimizerScales[1] =  lrTranslationScale;
  optimizerScales[2] =  lrTranslationScale;
  optimizerScales[3] =  lrTranslationScale;
  optimizerScales[4] =  lrTranslationScale;

  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximumStepLength( 1.0 ); 
  optimizer->SetMinimumStepLength( 0.0000001 );
  optimizer->SetNumberOfIterations( 100 );

  try { 
    lrRegistration->StartRegistration(); 
  } catch( itk::ExceptionObject & err )  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl;  
    return EXIT_FAILURE;
  } 

  OptimizerType::ParametersType lrParameters = 
          lrRegistration->GetLastTransformParameters();

  lrAngle           = lrParameters[0];
  lrRotationCenterX = lrParameters[1];
  lrRotationCenterY = lrParameters[2];
  lrTranslationX    = lrParameters[3];
  lrTranslationY    = lrParameters[4];

  // Print out results
  //
  std::cout << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "TRANSFORM COMPUTATION FOR LEFT/RIGHT ALIGN" << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " Angle (radians) " << lrAngle  << std::endl;
  std::cout << " Angle (degrees) " << lrAngle * 45.0 / atan(1.0)  << std::endl;
  std::cout << " Center X      = " << lrRotationCenterX  << std::endl;
  std::cout << " Center Y      = " << lrRotationCenterY  << std::endl;
  std::cout << " Translation X = " << lrTranslationX  << std::endl;
  std::cout << " Translation Y = " << lrTranslationY  << std::endl;
  std::cout << "==================================================" << std::endl;

/**
  // write images
  {
  typedef  unsigned short  OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< FixedImageType, OutputImageType > 
                    CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  filename = "leftBF.png";
  writer->SetFileName( filename );    
  caster->SetInput( leftPad->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  filename = "rightBF.png";
  writer->SetFileName( filename );    
  caster->SetInput( rightPad->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();
  }
**/

  // create transform to later apply to all images on the right hand side
  LRTransformType::Pointer lrTransform = LRTransformType::New();
  lrTransform->SetParameters( lrParameters );
 
  //=================================================================
  //=================================================================
  // Main loop registering frames
  //=================================================================
  //=================================================================

  // setup some stuff for the deformable registration we'll do between
  // frames, using the original rigid transform just to align the
  // left and right subframes
  
  const unsigned int SpaceDimension = Dimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;

  // specific transform and optimizer for deformable registration
  typedef itk::BSplineDeformableTransform< CoordinateRepType,
    SpaceDimension, SplineOrder >     DefTransformType;
  typedef itk::LBFGSOptimizer       DefOptimizerType;

  // make a registration object for deformable
  RegistrationType::Pointer defRegistration = RegistrationType::New();
  DefOptimizerType::Pointer   defOptimizer = DefOptimizerType::New();

  defRegistration->SetMetric       (metric);
  defRegistration->SetOptimizer    (defOptimizer);
  defRegistration->SetInterpolator (interpolator);

  DefTransformType::Pointer defTransform = DefTransformType::New();
  defRegistration->SetTransform(defTransform);

  // setting up the bsplinedeformabletransform grid.
  typedef DefTransformType::RegionType BSDTRegionType;
  BSDTRegionType bsplineRegion;
  BSDTRegionType::SizeType gridSizeOnImage;
  BSDTRegionType::SizeType gridBorderSize;
  BSDTRegionType::SizeType totalGridSize;

  gridSizeOnImage.Fill(20);
  gridBorderSize.Fill(3);
  totalGridSize = gridSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize(totalGridSize);

  typedef DefTransformType::SpacingType SpacingType;
  SpacingType spacing = leftWholeFixedExtractor->GetOutput()->GetSpacing();

  typedef DefTransformType::OriginType OriginType;
  OriginType origin = leftWholeFixedExtractor->GetOutput()->GetOrigin();

  FixedImageType::RegionType fixedRegion = 
    leftWholeFixedExtractor->GetOutput()->GetBufferedRegion();
  FixedImageType::SizeType fixedImageSize = fixedRegion.GetSize();

  for(unsigned int r=0; r<Dimension; r++)
    {
      spacing[r] *= floor( static_cast<double>(fixedImageSize[r] - 1)  / 
			   static_cast<double>(gridSizeOnImage[r] - 1) );
      origin[r]  -=  spacing[r]; 
    }
  
  defTransform->SetGridSpacing( spacing );
  defTransform->SetGridOrigin( origin );
  defTransform->SetGridRegion( bsplineRegion );
  
  typedef DefTransformType::ParametersType DefParametersType;
  const unsigned int defNumberOfParameters =
    defTransform->GetNumberOfParameters();

  DefParametersType defParameters(defNumberOfParameters);

  //
  // main loop
  //
  for (int framenum = startframe; framenum < endframe; framenum++) {    
    // create frame filename
    filename = buildFilename(inDir,fnameBase,framenum,fnameExt);
    
    cout << "Opening " << filename << endl;

    // set filename for moving image to current frame
    movingImageReader->SetFileName( filename );

    // preprocessing
    //    movingPreFilter->SetInput(movingImageReader->GetOutput());
    
    // set the inputs for the left/right side extractors
    /*
    leftMovingExtractor->SetInput(movingImageReader->GetOutput());
    leftMovingExtractor->Update();
    rightMovingExtractor->SetInput(movingImageReader->GetOutput());
    rightMovingExtractor->Update();
    */

    // set the inputs for the left/right whole side extractors
    leftWholeMovingExtractor->SetInput(movingImageReader->GetOutput());
    leftWholeMovingExtractor->Update();
    rightWholeMovingExtractor->SetInput(movingImageReader->GetOutput());
    rightWholeMovingExtractor->Update();

    
    // add files to registration object : we register on the left side
    defRegistration->SetFixedImage(  leftWholeFixedExtractor->GetOutput()  );
    defRegistration->SetMovingImage( leftWholeMovingExtractor->GetOutput() );
    fixedImageReader->Update();
    movingImageReader->Update();
    
    defRegistration->SetFixedImageRegion( 
            leftWholeFixedExtractor->GetOutput()->GetBufferedRegion() );

    // set up registration parameters

    defParameters.Fill(0.0);
    defTransform->SetParameters(defParameters);
    defRegistration->SetInitialTransformParameters(defTransform->GetParameters());
    
    std::cout << "Initial Deformable Parameters = " << std::endl;
    std::cout << defTransform->GetParameters() << std::endl;

    // setup optimizer
    defOptimizer->SetGradientConvergenceTolerance( 0.05 );
    defOptimizer->SetLineSearchAccuracy( 0.9 );
    defOptimizer->SetDefaultStepLength( 1.5 );
    defOptimizer->TraceOn();
    defOptimizer->SetMaximumNumberOfFunctionEvaluations( 1000 );
    
    itk::TimeProbesCollectorBase collector;
    
    std::cout << std::endl << "Starting Registration" << std::endl;

    try 
    { 
      collector.Start("Registration");
      defRegistration->StartRegistration(); 
      collector.Stop("Registration");
    } 
    catch( itk::ExceptionObject & err ) 
    { 
      std::cerr << "ExceptionObject caught during defReg !" << std::endl; 
      std::cerr << err << std::endl; 
      return EXIT_FAILURE;
    } 
    
    OptimizerType::ParametersType frameParameters = 
            defRegistration->GetLastTransformParameters();
    
    std::cout << "Last Transform Parameters" << std::endl;
    std::cout << frameParameters << std::endl;
    
    collector.Report();

    //=================================================================
    // OPTIMIZATION COMPLETE: EMIT REGISTERED IMAGE
    //=================================================================
    
    /*** 
    **** NOTE: WindowedSincInterpolateImageFunction still experimental!
    ***/
    /*
    const unsigned int WindowRadius = 5;
    
    // better interpolator for this purpose
    typedef itk::Function::HammingWindowFunction< WindowRadius > WindowFunctionType;

    // by default, the constant = 0 
    typedef itk::ConstantBoundaryCondition< MovingImageType > BoundaryConditionType;
    
    typedef itk::WindowedSincInterpolateImageFunction < MovingImageType, WindowRadius, 
        WindowFunctionType, BoundaryConditionType, double > 
      ResampleInterpolatorType;

    ResampleInterpolatorType::Pointer resampleInterpolator = 
      ResampleInterpolatorType::New(); 
    */
    
    // create resampling image filter
    typedef itk::ResampleImageFilter< MovingImageType, 
                                      FixedImageType > ResampleFilterType;
    
    // create transforms, one for each side of the image.  we need two
    // since the right hand side requires the lrTransform to be
    // composed in before it.
    LRTransformType::Pointer leftFrameTransform = LRTransformType::New();
    LRTransformType::Pointer rightFrameTransform = LRTransformType::New();

    defTransform->SetParameters(frameParameters);

// NO ADJUSTMENT FOR LR
    /*
    leftFrameTransform->SetParameters( frameParameters );
    rightFrameTransform->SetParameters( frameParameters );
    */

// APPLY LR ADJUSTMENT TO RIGHT

//    leftFrameTransform->SetParameters( frameParameters );
//    rightFrameTransform->SetParameters ( lrParameters );
    // compose lrTransform with rightFrame transform
    //    rightFrameTransform->Compose( defTransform, true);

// APPLY LR ADJUSTMENT TO LEFT
/*
    rightFrameTransform->SetParameters( frameParameters );
    leftFrameTransform->SetParameters ( lrTransform->GetParameters() );
    // compose lrTransform with rightFrame transform
    leftFrameTransform->Compose( regTransform, true);
*/
    ResampleFilterType::Pointer leftResampler = ResampleFilterType::New();
    ResampleFilterType::Pointer rightResampler = ResampleFilterType::New();
    
    //    leftResampler->SetTransform( leftFrameTransform );
    leftResampler->SetTransform( defTransform );
    leftResampler->SetInput( leftWholeMovingExtractor->GetOutput() );
    rightResampler->SetTransform( defTransform );
    rightResampler->SetInput( rightWholeMovingExtractor->GetOutput() );


    /**** UNCOMMENT BELOW TO GET DIFFERENT INTERPOLATOR ****/
    //    rightResampler->SetInterpolator( resampleInterpolator );
    //    leftResampler->SetInterpolator( resampleInterpolator );
    
    // left image
    FixedImageType::Pointer leftFixedImage = leftWholeFixedExtractor->GetOutput();
    
    leftResampler->SetSize( leftFixedImage->GetLargestPossibleRegion().GetSize() );
    leftResampler->SetOutputOrigin( leftFixedImage->GetOrigin() );
    leftResampler->SetOutputSpacing( leftFixedImage->GetSpacing() );
    leftResampler->SetDefaultPixelValue( 100 ); 

    // right image
    FixedImageType::Pointer rightFixedImage = rightWholeFixedExtractor->GetOutput();
    
    rightResampler->SetSize( rightFixedImage->GetLargestPossibleRegion().GetSize() );
    rightResampler->SetOutputOrigin( rightFixedImage->GetOrigin() );
    rightResampler->SetOutputSpacing( rightFixedImage->GetSpacing() );
    rightResampler->SetDefaultPixelValue( 100 ); 
    
    typedef  unsigned short  OutputPixelType;
    
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    
    typedef itk::CastImageFilter< FixedImageType, OutputImageType > 
                      CastFilterType;
    
    typedef itk::ImageFileWriter< OutputImageType >  WriterType;
    
    
    WriterType::Pointer      writer =  WriterType::New();
    CastFilterType::Pointer  caster =  CastFilterType::New();
    
    filename = buildFilename(outDir,"left"+fnameBase,framenum,fnameExt);
    writer->SetFileName( filename );    
    caster->SetInput( leftResampler->GetOutput() );
    writer->SetInput( caster->GetOutput()   );
    writer->Update();

    filename = buildFilename(outDir,"right"+fnameBase,framenum,fnameExt);
    writer->SetFileName( filename );    
    caster->SetInput( rightResampler->GetOutput() );
    writer->SetInput( caster->GetOutput()   );
    writer->Update();
    
  }
  
  return EXIT_SUCCESS;
}
