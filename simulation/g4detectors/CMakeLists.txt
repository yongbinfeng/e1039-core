# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(g4detectors CXX C)

# source code
include_directories("$ENV{OFFLINE_MAIN}/include/" "${PROJECT_SOURCE_DIR}/")
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

# ROOT
find_program(ROOTCONF "root-config")
if(ROOTCONF)
  message("-- Detecting ROOT:    found at ${ROOTCONF}")
else()
  message(FATAL_ERROR "-- Detecting ROOT:    not found")
endif()
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOT_PREFIX  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --evelibs   OUTPUT_VARIABLE ROOT_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Geant4
  find_program(G4CONF "geant4-config")
  if(G4CONF)
    message("-- Detecting Geant4:  found at ${G4CONF}")
  else()
    message(FATAL_ERROR "-- Detecting Geant4:  not found")
  endif()
  find_package(Geant4 REQUIRED gdml vis_all)
  execute_process(COMMAND geant4-config --prefix OUTPUT_VARIABLE GEANT4_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND geant4-config --cflags OUTPUT_VARIABLE GEANT4_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND geant4-config --libs   OUTPUT_VARIABLE GEANT4_LINK   OUTPUT_STRIP_TRAILING_WHITESPACE)

# ROOT dict generation
file(GLOB LinkDefhs ${PROJECT_SOURCE_DIR}/*LinkDef.h)
file(GLOB dicts "")
foreach(LinkDefh ${LinkDefhs})
	#message("LinkDefh: " ${LinkDefh})
	string(REPLACE "LinkDef.h" ".h" Dicth ${LinkDefh})
	string(REPLACE "${PROJECT_SOURCE_DIR}/" "" DictC ${LinkDefh})
	string(REPLACE "LinkDef.h" "_Dict.cc" DictC ${DictC})
	list(APPEND dicts ${DictC})
  #message("Dicth: " ${Dicth})
  #message("DictC: " ${DictC})
	add_custom_command(OUTPUT ${DictC} COMMAND rootcint ARGS -f ${DictC}
  -c
  -I$ENV{OFFLINE_MAIN}/include 
  -I${GEANT4_PREFIX}/include 
  ${Dicth} ${LinkDefh} DEPENDS ${LinkDefh})
endforeach(LinkDefh)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I$ENV{OFFLINE_MAIN}/include/ -I${ROOT_PREFIX}/include/ -I${GEANT4_PREFIX}/include/ ")

add_library(g4detectors_io SHARED
  ${PROJECT_SOURCE_DIR}/PHG4Parameters.cc
  PHG4Parameters_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4ParametersContainer.cc
  PHG4ParametersContainer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4ParameterContainerInterface.cc
  PHG4ParameterContainerInterface_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4ParameterInterface.cc
  PHG4ParameterInterface_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockCellGeom.cc
  PHG4BlockCellGeom_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockCellGeomContainer.cc
  PHG4BlockCellGeomContainer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockGeom.cc
  PHG4BlockGeom_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockGeomv1.cc
  PHG4BlockGeomv1_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockGeomContainer.cc
  PHG4BlockGeomContainer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4Cell.cc
  PHG4Cell_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4Cellv1.cc
  PHG4Cellv1_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4Cellv2.cc
  PHG4Cellv2_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CellContainer.cc
  PHG4CellContainer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CellDefs.cc
  PHG4CellDefs_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeom.cc
  PHG4CylinderGeom_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeomv1.cc
  PHG4CylinderGeomv1_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeomv2.cc
  PHG4CylinderGeomv2_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeomv3.cc
  PHG4CylinderGeomv3_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeomv4.cc
  PHG4CylinderGeomv4_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeom_Siladders.cc
  PHG4CylinderGeom_Siladders_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeom_MAPS.cc
  PHG4CylinderGeom_MAPS_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderGeomContainer.cc
  PHG4CylinderGeomContainer_Dict.cc
  PHG4CylinderCell_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellv1.cc
  PHG4CylinderCellv1_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCell_MAPS.cc
  PHG4CylinderCell_MAPS_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellv2.cc
  PHG4CylinderCellv2_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellv3.cc
  PHG4CylinderCellv3_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellContainer.cc
  PHG4CylinderCellContainer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellGeom.cc
  PHG4CylinderCellGeom_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellGeomContainer.cc
  PHG4CylinderCellGeomContainer_Dict.cc
)

add_library(g4detectors SHARED
  ${PROJECT_SOURCE_DIR}/PHG4BlockCellReco.cc
  PHG4BlockCellReco_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4BlockSubsystem.cc
  PHG4BlockSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/DPDigitizer.cc
  DPDigitizer_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4ConeDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4ConeRegionSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4ConeSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4ConeSubsystem.cc
  PHG4ConeSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderSubsystem.cc
  PHG4CylinderSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderCellReco.cc
  PHG4CylinderCellReco_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CylinderSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4PolyconDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4PolyconSubsystem.cc
  PHG4PolyconSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4PolyconSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4TargetCoilDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4TargetCoilSubsystem.cc
  PHG4TargetCoilSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4TargetCoilSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4BNLTargetCoilDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4BNLTargetCoilSubsystem.cc
  PHG4BNLTargetCoilSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4BNLTargetCoilSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4CollimatorDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4CollimatorSubsystem.cc
  PHG4CollimatorSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4CollimatorSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4E1039InsensDetector.cc
  ${PROJECT_SOURCE_DIR}/PHG4E1039InsensSubsystem.cc
  PHG4E1039InsensSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4E1039InsensSteppingAction.cc
  ${PROJECT_SOURCE_DIR}/PHG4DetectorGroupSubsystem.cc
  PHG4DetectorGroupSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4DetectorSubsystem.cc
  PHG4DetectorSubsystem_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4EventActionClearZeroEdep.cc
  ${PROJECT_SOURCE_DIR}/PHG4GenHit.cc
  PHG4GenHit_Dict.cc
  ${PROJECT_SOURCE_DIR}/PHG4StepStatusDecode.cc
  ${PROJECT_SOURCE_DIR}/PHG4TrackStatusDecode.cc
)

target_link_libraries(g4detectors_io -L./ -L$ENV{OFFLINE_MAIN}/lib -lphool -lboost_filesystem -lboost_system -lpdbcalBase -lphparameter -lg4testbench -L${ROOT_PREFIX}/lib -lXMLIO)
target_link_libraries(g4detectors    -L./ -lg4detectors_io -L$ENV{OFFLINE_MAIN}/lib -linterface_main -lgeom_svc -lphg4gdml -lphool -lSubsysReco -lg4testbench)

install(TARGETS g4detectors_io  					DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(TARGETS g4detectors     					DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

message(${CMAKE_PROJECT_NAME} " will be installed to " ${CMAKE_INSTALL_PREFIX})

file(GLOB dist_headers
  PHG4BlockCellGeom.h
  PHG4BlockCellGeomContainer.h
  PHG4BlockGeom.h
  PHG4BlockGeomContainer.h
  PHG4Cell.h
  PHG4Cellv1.h
  PHG4Cellv2.h
  PHG4CellContainer.h
  PHG4CellDefs.h
  PHG4CylinderCell.h
  PHG4CylinderCellv1.h
  PHG4CylinderCell_MAPS.h
  PHG4CylinderCellv2.h
  PHG4CylinderCellv3.h
  PHG4CylinderCellContainer.h
  PHG4CylinderGeom.h
  PHG4CylinderGeomv1.h
  PHG4CylinderGeomv2.h
  PHG4CylinderGeomv3.h
  PHG4CylinderGeomv4.h
  PHG4CylinderGeom_MAPS.h
  PHG4CylinderGeom_Siladders.h
  PHG4CylinderGeomContainer.h
  PHG4CylinderCellDefs.h
  PHG4CylinderCellGeom.h
  PHG4CylinderCellGeomContainer.h
  PHG4CylinderSubsystem.h
  PHG4PolyconSubsystem.h
  PHG4TargetCoilSubsystem.h
  PHG4BNLTargetCoilSubsystem.h
  PHG4CollimatorSubsystem.h
  PHG4E1039InsensSubsystem.h
  PHG4DetectorSubsystem.h
  PHG4DetectorGroupSubsystem.h
  PHG4ParameterContainerInterface.h
  PHG4ParameterInterface.h
  PHG4Parameters.h
  PHG4ParametersContainer.h
  PHG4StepStatusDecode.h
)

install(FILES ${dist_headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${CMAKE_PROJECT_NAME}/)


