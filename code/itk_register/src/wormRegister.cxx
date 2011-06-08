/***************************************************************/
/* wormRegister.cxx : An ITK-based tool for registering worm   */
/* movies from Shawn Lockery.                                  */
/*                                                             */
/* Matt Sottile / matt@cs.uoregon.edu                          */
/***************************************************************/

using namespace std;

#include <strings.h>
#include <sys/time.h>
#include <iostream>
#include "InputHandler.h"
#include "RigidWormRegister.h"
#include "sexp.h"

void usage(char *av0) {
  cout << "usage: " << av0 << " [config filename]" << endl;
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

void debug(string s) {
  cerr << "MAIN: " << s << endl;
}

/***********************************************************************
 * main()
 ***********************************************************************/
int main(int argc, char **argv) {
  //------------------------------------------------------------
  // parameter setup
  //------------------------------------------------------------
  
  if (argc != 2) usage(argv[0]);
    
  RegistrationParameters params(argv[1]);  
    
  //
  // setup registration
  //
  RigidWormRegistration *reg = new RigidWormRegistration();

  //
  // find current iteration of per-frame parameter log file
  //   (current_param_fname) and next iteration (new_param_file)
  //
  bool isLastParamfile = false;
  int i = 0;
  string current_param_fname;   
  string new_param_fname;   
  ifstream currentParamFile;
  
  while (!isLastParamfile)
  {  	 
    new_param_fname = buildFilename(params.output_dir, 
				    params.fname_base, i, "params");
    
    currentParamFile.open(new_param_fname.c_str(), ios::in);
    if(!currentParamFile.is_open())
      {
	if(i==0 && params.restarting) {
          cerr << "File '" << new_param_fname << "' not found. " <<
	    "If this is the first run, set the 'restarting' flag to " <<
	    "false in the config." << endl;
          exit(EXIT_FAILURE);
        }
	isLastParamfile = true;
	current_param_fname = buildFilename(params.output_dir, 
					    params.fname_base, --i, "params");
	
      }	
    else
      {
	i++;
	currentParamFile.close();
      }
  }	 
  
  //
  // extract optimization parameters from input file if they exist
  //
  RigidOptimizerParameters lrOptParams, ffOptParams;
  
  sexp_t *sx = params.getByName("left_right_params");
  if (sx != NULL)
    reg->extractParams(&lrOptParams, sx);
  else {
    cerr << "Using default left-right parameters." << endl;
    lrOptParams.scaleAngle = 1.0 / 100.0;
    lrOptParams.scaleRCX = lrOptParams.scaleRCY = 1.0 / 1000.0;
    lrOptParams.scaleTX = lrOptParams.scaleTY = 1.0 / 1000.0;
    lrOptParams.minStepLength = 0.0001;
    lrOptParams.maxStepLength = 1.0;
    lrOptParams.iterationCount = 200;
  }
  
  sx = params.getByName("frame_frame_params");
  if (sx != NULL)
    reg->extractParams(&ffOptParams, sx);
  else {
    cerr << "Using default frame-frame parameters." << endl;
    ffOptParams.scaleAngle = 1.0 / 25.0;
    ffOptParams.scaleRCX = ffOptParams.scaleRCY = 1.0 / 1000.0;
    ffOptParams.scaleTX = ffOptParams.scaleTY = 1.0 / 1000.0;
    ffOptParams.minStepLength = 0.0001;
    ffOptParams.maxStepLength = 1.0;
    ffOptParams.iterationCount = 150;
  }
  
  //
  // when in restart mode, take frame-to-frame parameters from the 
  //   log file of the last run 
  //  STEPLENGTH AND ITERATION COUNT ARE TAKEN FROM INPUT CONFIG
  //  
  if (params.restarting) {
    currentParamFile.open(current_param_fname.c_str(), ios::in);	  
    string str;
    int frameNum;
    while ( currentParamFile ) {  
      getline(currentParamFile,str);
      frameNum = atoi(str.substr(0, str.find(" [")).c_str());
      if(frameNum == params.restart_frame-1) {	  
	cerr << "Using params from frame " << frameNum << " => " << 
	  str.substr(str.find('[')) << endl;   
	float param[6];
	char delims[] = " ,[]";
	char *result = NULL;
	int i = 0;
	result = strtok( const_cast<char *>(str.c_str()), delims );
	while( result != NULL && i < 6) {		  
	  param[i] = atof(result);
	  i++;
	  result = strtok( NULL, delims );
	}
	ffOptParams.scaleAngle = param [1];
        ffOptParams.scaleRCX = param [2];
	ffOptParams.scaleRCY = param [3]; 
        ffOptParams.scaleTX = param [4];
	ffOptParams.scaleTY = param [5];
      }
    }
  }
  
  //
  // build filename for first file
  //
  string filename = buildFilename(params.input_dir, 
				  params.fname_base, 
				  params.ref_frame, 
				  params.extension);
  
  //
  // REGISTRATION SETUP
  //

  // 
  // setup filters if necessary
  //
  if (params.threshFilter) {
    reg->enableThreshFilter(params.loThresh, params.hiThresh);
  }

  //
  // set the fixed filename and moving filename to the same file since we're
  // going to start off doing the left/right registration
  //
  reg->setFixedFilename(filename);
  reg->setMovingFilename(filename);
    
  //
  // set the default extractors
  //
  reg->setDefaultExtractors();
  
  //
  // configure the registration to be the fixed image left vs right
  // 
  reg->setupLRRegistration();
  
  //
  // pass in our optimization parameters
  //
  reg->setOptimizationParams(lrOptParams);
  
  //
  // enable an observer
  //
  if (params.verbose == true)
    reg->enableObserver();
  
  //
  // perform the LR registration
  //
  cerr << "Performing left/right subframe registration..." << endl;
  reg->doRegistration();

  //
  // get results back and store away for later
  //
  RigidWormRegistration::TransformType::ParametersType lrTransformParams;
  
  lrTransformParams = reg->getLastRegistrationTransform();

  RigidWormRegistration::TransformType::ParametersType 
    lastFrameTransformParams;

  //
  // MAIN LOOP
  //  
  reg->setOptimizationParams(ffOptParams);
  reg->initializeByGeometry(params.initializeByGeometry);
  
  ofstream newParamFile;  
  
  int start_frame = params.start_frame;
  if(params.restarting)
    start_frame = params.restart_frame;
  
  timeval tv;
  int elapsedSec = 0;
  int tv_last;
  gettimeofday(&tv, NULL);
  tv_last = tv.tv_sec;
  for (int frame = start_frame; 
       frame < params.end_frame; 
       frame++) {
            
    //
    // current filename from frame
    //
    string moving_fname = buildFilename(params.input_dir, 
					params.fname_base, 
					frame, 
					params.extension);
    
    //
    // set moving image to be filename
    //
    reg->setMovingFilename(moving_fname);
    
    //
    // tell registration engine that we're doing frame/frame registration
    //
    reg->setupCrossFrameRegistration();
    
    //
    // do the registration
    //
    cerr << frame << "/" << params.end_frame << " working... ";
    reg->doRegistration();
    cerr << reg->getOptimizerValue() << endl;
    
    lastFrameTransformParams = reg->getLastRegistrationTransform();
    
    string left_outfname = buildFilename(params.output_dir, 
      "left_"+params.fname_base, frame, params.extension);

    string right_outfname = buildFilename(params.output_dir, 
      "right_"+params.fname_base, frame, params.extension);
      
    reg->writeTransformedMoving(left_outfname, right_outfname,
      lastFrameTransformParams, lrTransformParams);	  
	
    // 
    // append transform parameters for this frame to log file 
    //
    newParamFile.open(new_param_fname.c_str(), ios::app);
    newParamFile << frame << " " << lastFrameTransformParams << "\n";
    newParamFile.close();
    
    //
    // set starting parameters to last frames'.
    //
    reg->setTransformParameters(lastFrameTransformParams);  
    gettimeofday(&tv, NULL);
    elapsedSec += tv.tv_sec - tv_last;
    tv_last = tv.tv_sec;
    
    float estimatedRemaining = 
      (float)((params.end_frame - start_frame) - (frame - start_frame)) *
      ((float)elapsedSec / (float)((frame - start_frame)+1));
    
    cout << "Estimated remaining time (min.): " << 
      estimatedRemaining / 60 << endl;
  }
  
  return EXIT_SUCCESS;
}
