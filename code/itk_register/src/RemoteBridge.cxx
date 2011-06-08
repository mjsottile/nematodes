/**
 * Code related to binding the ITK registration code to an external
 * GUI via a socket-based protocol.  Not optimal, but isolates ITK's
 * annoying CMake-based build system from others such as QT or Cocoa.
 * Also allows us to avoid perturbing the registration codebase with
 * GUI-specific features.
 */

#include "sexp.h"

