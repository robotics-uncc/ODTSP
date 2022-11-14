// A. Wolek, 12-May-2017

#ifndef ODTSP_TARGET_H
#define ODTSP_TARGET_H

#include<vector>
#include "ODTSP_Definitions.h"
#include "ODTSP_Problem.h"


namespace ODTSP {


enum TargetType {
  Single, // contains a single target
  Cluster // contains multiple targets
};

class Target {

private: 
  // constructor arguments:
  vi elems_;  // integer array of elements in the target 
              // the elements refer to targets listed in ODTSP::Problem
              // (a cluster has multiple elements)

  double centerX_, centerY_; // orbit center x,y 
  double orbitRadius_, orbitAngleRad_; // orbit radius and angle
  ODTSP::TargetType type_; // type of target (see above)

  // The motion along the orbit can be specified several ways: 
  // Option A: using set_entryState()
  double entryX_, entryY_, entryHrad_; // entryHrad measured (+) CCW from x-axis
  // Option B: using set_exitState()
  double exitX_, exitY_, exitHrad_;
  // Option C: set_entryAngleOrbitDir()
  double entryAngleRad_; // entry/exit polar angles measured (+) CCW from x-axis
  int orbitDir_; // (+1 = CCW, -1 = CW)
  // Regardless of the set method used, all of the above variables are calculated, including:
  double exitAngleRad_;

  // Useful conversion function 
  double polarAngleFromState(const double & x, const double & y);

  bool targetDefined_; // flag checking if the appropriate set method has been 
                      // involed to fully defined the object

  // error checking
  double tolAngleRad_; // tolerance for declaring angles are equal
  double tolLength_; // tolerance for declaring lengths are equal

  // check if a point is on the orbit (based on range to center)
  bool pointOnOrbit(const double & x,const double & y);

  // determine the orbit direction given a heading and polar angle,
  // also check if it is valid
  bool tangentHeading(const double & headingRad, 
                                   const double & polarAngleRad, 
                                   int & orbitDir);
  // check if angle is between [0,2pi)
  bool validAngle(const double & h);

public:
  // constructor, destructor
  Target(){};
  ~Target(){};

  // required to define a target 
  void initialize(const std::vector<int> & elems, const double & centerX, 
                  const double & centerY, const double & orbitRadius,
                  const double & orbitAngleRad, 
                  const ODTSP::TargetType & type);

  // optional modification to error checking defaults
  bool set_tolerance(const double & tolAngleRad, const double & tolLength);

  // the user must set one of these functions to fully define the object:
  // successive calls to these functions will over-write previous values
  bool set_entryState(const double & entryX, const double & entryY, 
                      const double & entryHrad);
  bool set_exitState(const double & exitX, const double & exitY, 
                     const double & exitHrad);
  bool set_entryAngleOrbitDir(const double & entryAngleRad, 
                              const int & orbitDir);

  // display object contents to screen
  void print(); 

  // the target must be defined for these functions to provide an output,
  // otherwise an exception is thrown 
  double get_entryX();
  double get_entryY();
  double get_entryHrad();
  double get_entryAngleRad();

  double get_exitX();
  double get_exitY();
  double get_exitHrad();
  double get_exitAngleRad();

  double get_centerX();
  double get_centerY();
  double get_orbitRadius();
  double get_orbitAngleRad();
  int get_orbitDir();
  vi get_elems();
  int get_numElems();
  bool contains_elem(int elemInd);
  bool contains_elem(int elemInd, int & entryInd);
  ODTSP::TargetType get_type();


};

} // namespace

#endif

