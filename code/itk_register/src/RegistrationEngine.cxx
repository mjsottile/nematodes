/******************************************************************/
/* RegistrationEngine.cxx                                         */
/******************************************************************/
#include "RegistrationEngine.h"

/*
 * constructor : instantiate fields and set flags to false.
 */
RegistrationEngine::RegistrationEngine() {
  // Set up image readers - one for fixed, one for moving frame
  fixedImageReader = ImageReaderType::New();
  movingImageReader = ImageReaderType::New();

  // Create two extractors per frame : one for the left half and one
  // for the right half.
  leftFixedExtractor = ExtractorType::New();
  leftMovingExtractor = ExtractorType::New();
  rightFixedExtractor = ExtractorType::New();
  rightMovingExtractor = ExtractorType::New();

  // Make an image writer
  imageWriter = ImageWriterType::New();
  
  // Flags
  fixedImageSet = false;
  movingImageSet = false;
  writerIsSet = false;
}

/*
 * destructor : does nothing for now.
 */
RegistrationEngine::~RegistrationEngine() {
  // empty  
}



void
RegistrationEngine::setupCrossFrameRegistration() {
  debug("DO NOT CALL REGISTRATION ENGINE DIRECTLY!");
  exit(EXIT_FAILURE);
}

void 
RegistrationEngine::doRegistration() {
  debug("DO NOT CALL REGISTRATION ENGINE DIRECTLY!");
  exit(EXIT_FAILURE);  
}

void
RegistrationEngine::setupLRRegistration() {
  debug("DO NOT CALL REGISTRATION ENGINE DIRECTLY!");
  exit(EXIT_FAILURE);
}



void 
RegistrationEngine::configureExtractors(ImageReaderType::Pointer reader,
					ExtractorType::Pointer leftExtractor,
					ExtractorType::Pointer rightExtractor) {
  // Force an update to read an image
  reader->Update();
            
  // Get the size of the image
  ImageType::SizeType imgSize =
    reader->GetOutput()->GetLargestPossibleRegion().GetSize();

  // Variables for holding shape and location info for the halves.
  ImageType::IndexType leftStart, rightStart;
  ImageType::SizeType leftSize, rightSize;

  // Left side
  leftStart[0] = 0;            leftStart[1] = 0;
  leftSize[0]  = imgSize[0]/2; leftSize[1]  = imgSize[1];

  // Right side
  rightStart[0] = leftSize[0];   rightStart[1] = 0;
  rightXOffset  = rightStart[0]; rightYOffset  = rightStart[1];
  rightSize[0]  = imgSize[0]/2;  rightSize[1]  = imgSize[1];

  // Create regions out of shape/location info.
  ImageType::RegionType leftRegion, rightRegion;

  leftRegion.SetSize(leftSize);
  leftRegion.SetIndex(leftStart);
  rightRegion.SetSize(rightSize);
  rightRegion.SetIndex(rightStart);

  // Use these regions as the ROIs for each extractor
  leftExtractor->SetRegionOfInterest(leftRegion);
  rightExtractor->SetRegionOfInterest(rightRegion);
}

void RegistrationEngine::setFixedFilename(string fname) {
  fixedImageReader->SetFileName(fname);
  fixedImageReader->Update();
  fixedImageSet = true;

  configureExtractors(fixedImageReader, leftFixedExtractor, 
		      rightFixedExtractor);
}

void RegistrationEngine::setMovingFilename(string fname) {
  movingImageReader->SetFileName(fname);
  movingImageReader->Update();
  movingImageSet = true;

  configureExtractors(movingImageReader, leftMovingExtractor,
		      rightMovingExtractor);
		      
}

void RegistrationEngine::writeImage(string fname, SourceType::Pointer src) {
  CastFilterType::Pointer caster = CastFilterType::New();
  
  imageWriter->SetFileName(fname);
  caster->SetInput(src->GetOutput());
  imageWriter->SetInput(caster->GetOutput());
  imageWriter->Update();
}
