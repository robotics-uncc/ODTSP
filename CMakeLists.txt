project (ODTSP CXX C)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND cmake_policy)
	CMAKE_POLICY(SET CMP0003 OLD)
	CMAKE_POLICY(SET CMP0005 OLD)
ENDIF(COMMAND cmake_policy)

find_package(Ceres REQUIRED) # Note: Ceres must have been installed with make install.
include_directories(${CERES_INCLUDE_DIRS})


# DTSP (Jim)
ADD_DEFINITIONS(-DARCH_K8)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ADD_DEFINITIONS(-Wunused-result)

# set flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC -g") #-std=c++11 -fPIC")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# get path variables
get_filename_component(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

# include headers into the project
include_directories(
# external
${ROOT_DIR}/ODTSP/external/miniball/
# custom packages
${ROOT_DIR}/ODTSP/src
)

# external libraries
link_directories(
/lib
/usr/lib
/usr/local/lib
)

#############################################################################

# use function file(GLOB...) to grab all source files using wildcard syntax
file(GLOB SOURCES_BASIC "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
set(SOURCES ${SOURCES_BASIC})

#set(LIBS "RobustDubins;MathTools;ConvectedDubins")
set(LIBS ${CERES_LIBRARIES})

message(STATUS "SOURCES: ${SOURCES}")
message(STATUS "LIBS: ${LIBS}")

add_library(ODTSP ${SOURCES})
target_link_libraries(ODTSP "${LIBS}")

# build programs
file(GLOB PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/programs/*.cpp")  
foreach( testsourcefile ${PROGRAMS} )
  # I used a simple string replace, to cut off .cpp.
  get_filename_component(testname ${testsourcefile} NAME_WE)
  add_executable(${testname} ${testsourcefile} )
  set_target_properties(${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "../bin")
  # Make sure YourLib is linked to each app
  target_link_libraries(${testname} ODTSP)
endforeach()

# build programs
file(GLOB PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/programs/dopt/*.cpp")  
foreach( testsourcefile ${PROGRAMS} )
  # I used a simple string replace, to cut off .cpp.
  get_filename_component(testname ${testsourcefile} NAME_WE)
  add_executable(${testname} ${testsourcefile} )
  set_target_properties(${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "../bin")
  # Make sure YourLib is linked to each app
  target_link_libraries(${testname} ODTSP)
endforeach()



