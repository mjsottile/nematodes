cmake_minimum_required(VERSION 2.6)
PROJECT(wormRegister)

FIND_PACKAGE(ITK)
IF (ITK_FOUND)
  INCLUDE(${USE_ITK_FILE})
ELSE (ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK Not Found!  Please set ITK_DIR.")
ENDIF (ITK_FOUND)

INCLUDE_DIRECTORIES(
  /Users/matt/Research/sexpr/sexpr/src
)

INCLUDE_DIRECTORIES(
  ${wormRegister_SOURCE_DIR}
)

LINK_DIRECTORIES(
  /Users/matt/Research/sexpr/sexpr/src
)

ADD_EXECUTABLE( wormRegister 
  wormRegister.cxx 
  RigidWormRegister.cxx 
  FilterOperations.cxx
  RegistrationEngine.cxx
  InputHandler.cxx
)

TARGET_LINK_LIBRARIES ( wormRegister
  ITKBasicFilters 
  ITKCommon 
  ITKIO
  ITKFEM
  ITKNumerics
  sexp
)
