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

#ifndef CONVECTED_DUBINS_PATH_H
#define CONVECTED_DUBINS_PATH_H

#include<vector>
#include<string>

#include "ConvectedDubins_Problem.h"

namespace ConvectedDubins {


// The ConvectedDubins::Path object can be used on its own (e.g., to draw a 
// path with known parameters, or it can be the output from 
// the ConvectedDubins::Solver

enum pathClass { BSB, BBB, NO_CLASS };

enum pathType { LSL, RSR, LSR, RSL, LRL, RLR, NO_TYPE };

enum pathStatus { INFEASIBLE, FEASIBLE, OPTIMAL, NO_STATUS };

class Path {

private: 

  // Note: BSB paths are defined with
  // 1st segment: trochoid with t \in [0, tA)
  // 2nd segment: straight line with t \in [tA, tBeta)
  // 3rd segment: trochoid with t\in [tBeta, T]
  pathClass pc_;
  pathType pt_;
  pathStatus ps_;   

  // to avoid storing redundant information, the ConvectedDubins::Path can 
  // reference an already defined problem 
  ConvectedDubins::Problem * prob_;
  bool probSuppliedFlag; // flag indicating if prob_ pointer supplied 
  ConvectedDubins::Problem probGen_; // used only when a prob pointer is not 
                                     // given. In this case the user 
                                     // provides initial state, final state, etc.
                                     // via this Path object rather than via 
                                     // the Problem object 

  // These parameters define a Convected Dubins path 
  // segment A: trochoid 
  double tA_; // time duration
  double dirA_; // direction 

  // segment B: straight, or trochoid 
  double tBeta_; // elapsed time to start trochoid C  

  // segment C: trochoid 
  double T_; // elapsed total time 
  double dirB_; // direction 
 
  // plotting
  double dt_; // time-step for defining state history
  vd x_, y_, h_, t_;  // state history in the inertial frame
  bool pathHistoryComputedFlag_; // flag indicating computePathHistory() was run 
 
  // constants 
  double t2pi_;

  // private set functions
  void initialize();

  // intermediate states at junction of extremals
	double xFinal_; // Final x position, default=0 (LU)
	double yFinal_; // Final y position, default=0 (LU)
	double hFinalRad_; // Final heading angle, default=0 (rad)
  double xBinit_, yBinit_, hBinitRad_; // at start of B extremal 
  double xBfinal_, yBfinal_, hBfinalRad_; // at end of B extremal 
  bool connectedStatesComputedFlag_; // flag indicating if these intermediate
                                     // state have been computed 

  // function returns the state at a particular time 
  void pathXYH( const double & t, double & x, double & y, double & h);

public: 
  // constructor and destruct
	Path();
  Path( ConvectedDubins::Problem * prob );
	virtual ~Path(){};

  // if * prob was not specified in constructor, the following must be called:
  void set_prob( ConvectedDubins::Problem * prob  );
  void set_windMagXY( double windMagX, double windMagY);
  void set_windMagDir( double windMagnitude, double windDirRad );
  void set_vehicleProperties( double minTurningRadius, double nominalSpeed );
  void set_stateInitial( const vd & stateInitial );
  void set_stateInitial( double xInitial, double yInitial, double hInitialRad );

  // once the initial state and other problem parameters are defined (as above) 
  // the path can be specified using its parametrization (directions, and 
  // durations of each extremal segment)

  // paths defined in this way are automatically elevated to FEASIBLE status
  void set_pathParamsBSB( double dirA, double dirB, double tA, double tBeta, double T);
  void set_pathParamsBBB( double dirA, double tA, double tBeta, double T);

  // a path is set optimal, typically by the solver 
  void set_pathStatusOptimal();

  // similarily, a path is set infeasible, typically by the solver
  void set_pathStatusInfeasible();

  // set the type directly
  void set_pathType( ConvectedDubins::pathType pt){ pt_ = pt;};

  // print data to screen
  void print();  

  // path must be defined if it is FEASIBLE or OPTIMAL
  bool isDefined();

  // use 'del' displacment functions to compute the states at junctions
  // of extremals 
  void computeConnectingStates();

  // the path history is computed by discretized the elapsed time interval 
  // [0, T] and recursively calling the pathXYH function
  // the result is stored in the x_, y_ and h_ vector variables
  void computePathHistory();

  // print wpts to screen 
  void printWaypoints();

  // write X,Y coordinates of path to a file compatible with Octave
  // (and MATLAB with small changes)
  void plotXY(const std::string & fileName, const int & figureNumber);

  // get functions 
  ConvectedDubins::pathStatus get_pathStatus(){return ps_;};
  ConvectedDubins::pathClass get_pathClass(){return pc_;};
  ConvectedDubins::pathType get_pathTyped(){return pt_;};
  double get_T(){return T_;};
  double get_xFinal(){return xFinal_;};
  double get_yFinal(){return yFinal_;};
  double get_hFinalRad(){return hFinalRad_;};
  void get_waypoints( vd & x, vd & y, vd & hRad );
  void get_waypoints_XYEastNorth( vd & x, vd & y, vd & hRad );

}; // class

} // namespace
#endif
