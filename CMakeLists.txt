cmake_minimum_required(VERSION 2.8)
#-----------------------------------------------------------------------------
set(MODULE_NAME QIsModule)

find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})

find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

set(QISMODULE_SOURCE QIs.cxx)
GENERATECLP(QISMODULE_SOURCE QIs.xml)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/ITKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
