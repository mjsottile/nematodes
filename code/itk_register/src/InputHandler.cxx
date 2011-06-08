#include "InputHandler.h"
#include <iostream>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

RegistrationParameters::RegistrationParameters(string fname) {
  int fd;
  sexp_iowrap_t *iow;
  
  // set defaults
  input_dir = ".";
  output_dir = "/tmp/";
  extension = "png";
  fname_base = "";
  fnum_digits = 4;
  ref_frame = 1;
  start_frame = 2;
  restart_frame = 2;
  end_frame = 3;
  restarting = false;
  verbose = false;
  initializeByGeometry = false;
  threshFilter = false;
  loThresh = 0;
  hiThresh = 100000;
  
  fd = open(fname.c_str(), O_RDONLY);
  
  iow = init_iowrap(fd);
  
  param_sx = read_one_sexp(iow);
  
  if (param_sx == NULL) {
    cerr << "ERROR READING PARAMETERS!" << endl;
    exit(EXIT_FAILURE);
  }
  
  extract();
  
  destroy_iowrap(iow);
  
  close(fd);
}
  
RegistrationParameters::~RegistrationParameters() {
  if (param_sx != NULL) {
    destroy_sexp(param_sx);
    param_sx = NULL;
  }
}

string RegistrationParameters::getOneParam(string pname, string def) {
  sexp_t *param;
  char pstr[1024];
  string s;
  
  param = find_sexp(pname.c_str(), param_sx);
  
  if (param == NULL) {
    return def;
  }
  
  if (param->next->ty == SEXP_VALUE) {
    strncpy(pstr,param->next->val,param->next->val_used);
    pstr[param->next->val_used] = '\0';
  } else {
    print_sexp(pstr,1024,param->next);
  }
  
  s = pstr;
  return s;
}
  
sexp_t 
*RegistrationParameters::getByName(string name) 
{
  sexp_t *sx = find_sexp(name.c_str(), param_sx);
  return sx;
}
  
void RegistrationParameters::extract() 
{
  string tmp;
  
  input_dir = getOneParam("input_dir", input_dir);
  output_dir = getOneParam("output_dir", output_dir);
  extension = getOneParam("extension", extension);
  fname_base = getOneParam("filename_base", fname_base);
  tmp = getOneParam("frame_digits", "");
  if (tmp != "") {
    fnum_digits = atoi(tmp.c_str());
  }
  tmp = getOneParam("reference_frame", "");
  if (tmp != "") {
    ref_frame = atoi(tmp.c_str());
  }
  tmp = getOneParam("start_frame", "");
  if (tmp != "") {
    start_frame = atoi(tmp.c_str());
  }
  tmp = getOneParam("restart_frame", "");
  if (tmp != "") {
    restart_frame = atoi(tmp.c_str());
  }
  tmp = getOneParam("end_frame", "");
  if (tmp != "") {
    end_frame = atoi(tmp.c_str());
  }
  tmp = getOneParam("restarting", "false");
  if (tmp == "true") {
    restarting = true;
  }
  tmp = getOneParam("verbose", "false");
  if (tmp == "true") {
    verbose = true;
  }
  tmp = getOneParam("initializer", "mass");
  if (tmp == "geometry") {
    initializeByGeometry = true;
  }
  tmp = getOneParam("threshold_filter", "false");
  if (tmp == "true") {
    threshFilter = true;
    tmp = getOneParam("lowerThreshold", "");
    if (tmp != "") {
      loThresh = atoi(tmp.c_str());
    }
    tmp = getOneParam("upperThreshold", "");
    if (tmp != "") {
      hiThresh = atoi(tmp.c_str());
    }
  }
  
  if (verbose) {
    cout << "=====[ BASE PARAMETER SET ]=====" << endl;
    cout << "INPUT DIRECTORY: " << input_dir << endl;
    cout << "OUTPUT DIRECTORY: " << output_dir << endl;
    cout << "EXTENSION: " << extension << endl;
    cout << "FILENAME BASE: " << fname_base << endl;
    cout << "FRAME DIGITS: " << fnum_digits << endl;
    cout << "REF/START/RESTART/END: " << ref_frame << "/" << start_frame << "/" <<
      restart_frame << "/" << end_frame << endl;
    cout << "INITIALIZE BY GEOMETRY: " << initializeByGeometry << endl;
	cout << "RESTARTING: " << restarting << endl;
    cout << "VERBOSE: " << verbose << endl;
    cout << "THRESHOLD FILTER: " << threshFilter;
    if (threshFilter) {
      cout << " LO:" << loThresh << "  HI: " << hiThresh << endl;
    } else {
      cout << endl;
    }
  }
}
