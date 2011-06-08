/* ======================================================================

   RegistrationEngine class for worm registration code.  This is
   intended to be the base class for both rigid and deformable 
   techniques, encapsulating common operations like reading, writing,
   and extracting regions of interest from the images.
   
   Matthew Sottile / matt@cs.uoregon.edu
   Summer 2008

   ====================================================================== */
#ifndef __REGISTRATIONENGINE_H__
#define __REGISTRATIONENGINE_H__

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"


const unsigned int Dimension = 2;

using namespace std;

/**
 * RegistrationEngine class.
 *
 * This is a virtual base class that contains the basics that all
 * registration pipelines will need, such as the image reader and
 * writer, basic PixelType and ImageTypes, and ROI Filters for
 * extracting specific parts of an image for registration.  This class
 * specifically assumes that each frame (both moving and fixed) may have
 * two meaningful halves that should be registered separately, hence the
 * left/right separation of extractors.
 */
class RegistrationEngine
{
 //
 // typedefs
 //
 public:
  // assume 16-bit pixel types
  typedef unsigned short                      PixelType;
  typedef unsigned short                      OutputPixelType;

  typedef itk::Image<PixelType, Dimension>       ImageType;
  typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
  
  typedef itk::ImageFileReader<ImageType>          ImageReaderType;
  typedef itk::ImageFileWriter<OutputImageType>    ImageWriterType;

  typedef itk::RegionOfInterestImageFilter<ImageType, ImageType>
    ExtractorType;

  typedef itk::ImageSource<ImageType> SourceType;
  
  typedef itk::ResampleImageFilter<ImageType, ImageType>
    ResampleFilterType;

  typedef itk::CastImageFilter<ImageType, OutputImageType> CastFilterType;

 //
 // public member functions
 //
 public:

  RegistrationEngine();
  virtual ~RegistrationEngine();

  void setFixedFilename(string fname);
  void setMovingFilename(string fname);

  virtual void doRegistration();

  void configureExtractors(ImageReaderType::Pointer reader,
			   ExtractorType::Pointer leftExtractor,
			   ExtractorType::Pointer rightExtractor);
			   
  virtual void setupLRRegistration();
  virtual void setupCrossFrameRegistration();     

  void writeImage(string fname, SourceType::Pointer src);

 //
 // protected fields that the public shouldn't see but children
 // who inherit from this class should.
 //
 protected:

  ImageReaderType::Pointer     fixedImageReader;
  ImageReaderType::Pointer     movingImageReader;

  bool fixedImageSet, movingImageSet;
  bool writerIsSet;

  int rightXOffset, rightYOffset;

  ExtractorType::Pointer       leftFixedExtractor;
  ExtractorType::Pointer       leftMovingExtractor;
  ExtractorType::Pointer       rightFixedExtractor;
  ExtractorType::Pointer       rightMovingExtractor;

  ImageWriterType::Pointer     imageWriter;
  
private:
  void debug(string s) { cerr << "RegistrationEngine: " << s << endl;}

};

#endif /* __REGISTRATIONENGINE_H__ */
