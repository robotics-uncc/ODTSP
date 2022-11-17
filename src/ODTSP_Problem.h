// Orbiting Dubins Traveling Salesman Problem
// A. Wolek, J. McMahon, Sep-2018
//
// LICENSE:
//
// The source code is in the public domain and not licensed or under
// copyright. The information and software may be used freely by the public.
// As required by 17 U.S.C. 403, third parties producing copyrighted works
// consisting predominantly of the material produced by U.S. government
// agencies must provide notice with such work(s) identifying the U.S.
// Government material incorporated and stating that such material is not
// subject to copyright protection.
//
// Derived works shall not identify themselves in a manner that implies an
// endorsement by or an affiliation with the Naval Research Laboratory.
//
// RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
// SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
// RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
// OF RECIPIENT IN THE USE OF THE SOFTWARE.
//
// (UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

#ifndef ODTSP_PROBLEM
#define ODTSP_PROBLEM

#include "ODTSP_Definitions.h"
#include "RobustDubins_Problem.h"
#include "RobustDubins_Solver.h"
#include "RobustDubins_Path.h"
#include "MathTools.h"

namespace ODTSP {

enum motionModel { DUBINS, CONVECTED_DUBINS, FEASIBLE_DUBINS, UNDEFINED_MOTION_MODEL };

class Problem {

private:
  // required inputs
  int numTargs_; // the number of targets to revisit
  vd xTargs_; // x location of the targets
  vd yTargs_; // y location of the targets
  double orbitAngleRad_; // angle subtended by the orbit
  double orbitRadius_; // minimum distance for effective sonar use 
                       // (radius of the RID circle for a single target)
  double orbitRadiusMax_; // maximum distance for effective sonar use
  double aspectAngleMaxRad_; // max angle for viewing target (rel. broadside) 

  // vehicle parameters  
  ODTSP::motionModel mm_; // enumeration, see above 
	double minTurnRadius_; // minimum turning radius R 
  double v_; // vehicle speed

  // convected dubins
  double wx_, wy_; // wind speed in x and y direction
  double w_, hwRad_; // wind magnitude and direction in ENU frame
  double hwRadNED_; // wind magnitude converted to NED Frame 

  // derived 
  double dmax_; // computed reference value

public: 
  // define a Dubins problem 
  Problem();

  // define a Convected Dubins problem 
  //Problem();

  ~Problem(){};

  // required
  bool set_targs( const vd & xTargs, const vd & yTargs );

  // 
  bool set_orbitProperties( double orbitAngleRad, double orbitRadius );

  // optional (for clustering)
  bool set_orbitProperties( double orbitAngleRad, double orbitRadius,
                            double orbitRadiusMax, double aspectAngleMaxRad );

  // use if motion model is DUBINS
  bool set_motionProperties( double minTurnRadius );

  // use if motion model is CONVECTED_DUBINS
  bool set_motionProperties( double minTurnRadius, double vehicleSpeed ,  ODTSP::motionModel mm );
  
  // use one of the following if motione model is CONVECTED_DUBINS
  void set_windMagDir( double windMagnitude, double windDirRad);
  void set_windMagXY( double windMagX,  double windMagY);

  // main
  void print(); // prints the problem statement
  bool isDefined();

  // get functions
  int     get_numTargs(){return numTargs_;};
  vd      get_xTargs(){return xTargs_;};
  double  get_xTargs(int index){return xTargs_.at(index);};
  vd      get_yTargs(){return yTargs_;};
  double  get_yTargs(int index){return yTargs_.at(index);};
  double  get_orbitAngleRad(){return orbitAngleRad_;};
  double  get_orbitRadius(){return orbitRadius_;};

  // these must have been set with clustering option 'motionProperties'
  double  get_orbitRadiusMax();
  double  get_aspectAngleMaxRad();
  double  get_dmax();

  // get functions 
  double  R(){return minTurnRadius_;};
  double  v(){return v_;}; 
  double  wx(){return wx_;};
  double  wy(){return wy_;};
  double  w(){return w_;};
  double  hwRad(){return hwRad_;};
  double  hwRadNED(){return hwRadNED_;};

  ODTSP::motionModel motionModel(){return mm_;};


};

} // namespace

#endif

