#ifndef __INPUTHANDLER_H__
#define __INPUTHANDLER_H__

#include <string>
#include "sexp.h"

class RegistrationParameters {
public:
  std::string input_dir;
  std::string output_dir;
  std::string extension;
  std::string fname_base;
  int    fnum_digits;
  sexp_t *param_sx;
  int    ref_frame, start_frame, restart_frame, end_frame;
  bool   restarting;
  bool   verbose;
  bool   initializeByGeometry;
  bool   threshFilter;
  long   loThresh, hiThresh;
  
  RegistrationParameters(std::string fname);
  ~RegistrationParameters();
  void extract();
  sexp_t *getByName(std::string name);
  
private:
  std::string getOneParam(std::string pname, std::string def);
};

#endif /* __INPUTHANDLER_H__ */
