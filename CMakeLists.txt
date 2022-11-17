# Orbiting Dubins Traveling Salesman Problem
# A. Wolek, J. McMahon, Sep-2018
#
# LICENSE:
#
# The source code is in the public domain and not licensed or under
# copyright. The information and software may be used freely by the public.
# As required by 17 U.S.C. 403, third parties producing copyrighted works
# consisting predominantly of the material produced by U.S. government
# agencies must provide notice with such work(s) identifying the U.S.
# Government material incorporated and stating that such material is not
# subject to copyright protection.
#
# Derived works shall not identify themselves in a manner that implies an
# endorsement by or an affiliation with the Naval Research Laboratory.
#
# RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
# SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
# RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
# OF RECIPIENT IN THE USE OF THE SOFTWARE.
#
# (UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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



