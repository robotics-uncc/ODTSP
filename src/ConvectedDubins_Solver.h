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

#ifndef CONVECTED_DUBINS_SOLVER_H
#define CONVECTED_DUBINS_SOLVER_H

#include "ConvectedDubins_Problem.h"
#include "ConvectedDubins_Path.h"

namespace ConvectedDubins {

class Solver {

private: 
  ConvectedDubins::Problem * prob_;

  // tolerances
  int N_ = 20; // used to define length of "pvec" that stores the parameters
               // defining a solution. This length is the number of times 
               // the fixed point iterator runs. 
  
  double tol_ = 0.001; // defines a stopping criteria
                       // to exit before N_ iterations, if consecutive pvec
                       // values are within this tolerance 

  double eps_ = 0.001; // used for checking sign conditions with some tolerance
                      // eg. (x < -eps) checks if x is negative but allows 
                      // some tolerance for numerical errors 
                      // (e.g., if x = -0.0001)
                      // also used to check that x and y are equal

  int numTestPts_ = 10; // the fixed point solver is started 
                        // with this many start points spaced uniformly

  double sameRootEpsilon_ = 0.01; // most of the time the fixed point solver  
                                  // returns the same root from nearby start 
                                  // points. this values is used to check if 
                                  // two roots are the same

  double straightLineCourseTolRad_ = 0.01; // when the initial course angle 
                                           // points at the destination
                                           // and the initial and final heading
                                           // are equal than the solution is a 
                                           // straight line path. Here this 
                                           // tolerance is used to enter 
                                           // a special case where straight 
                                           // solution is a straight line.

  // initial and final state in the trochoidal frame 
  double x0_, y0_;
  double x1_, y1_;

  // intermediate variables used by solver
  double t2pi_; 

  // paths 
  ConvectedDubins::Path LSL_;
  ConvectedDubins::Path RSR_;
  ConvectedDubins::Path LSR_;
  ConvectedDubins::Path RSL_;
  ConvectedDubins::Path RLR_;
  ConvectedDubins::Path LRL_;

  // pointer to the optimal path   
  ConvectedDubins::Path * optPath_;

  // invokes either the analytical solution, or the fixedPointBSB
  void solveBSB( ConvectedDubins::pathType pt );

  // invokes the fixedPointBBB
  void solveBBB( ConvectedDubins::pathType pt );

  // solvers via contraction mapping 
  double fixedPointBSB( double p0, double d1, double d2, double k,
                        double xt10, double yt10, double phit1,
                        double xt20, double yt20, double phit2 );
  bool fixedPointBBB( std::vector<double> pvecInit, double d1, double k, 
                      std::vector<double> & pvecOut,
                      double xt10, double yt10, double phit1 );
  
  // check/set if a straight line path is the optimal solution 
  bool straightLineSolnExists();

  // verify that a candidate path is valid by checking several criteria 
  bool checkConditionsBSB( double d1, double d2, double tA, double tB, 
                           double xt10, double yt10, double phit1,
                           double xt20, double yt20, double phit2,
                           ConvectedDubins::pathType pt  );
  bool checkConditionsBBB( double d1, double tA, double tB, double T,
                           double xt10, double yt10, double phit1,
                           ConvectedDubins::pathType pt );

  // retrieves cost (elapsed time) for each path type, computes 
  // minimum and sets corresponding optimal 
  void setOptimal();

  
public: 

  // constructor and destruct
	Solver( ConvectedDubins::Problem * prob );
	virtual ~Solver(){};

  // run computeConnectingStates() functions for each path type to obtain
  // final endpoints, intermediate extremal junction points 
  void computeConnectingStates();

  // print solutions to screen 
  void printSolns();

  // write XY path coordinates of each path to Octave file 
  void plotSolns(std::string fileName);

  // get functions   
  ConvectedDubins::Path get_path( ConvectedDubins::pathType type );
  ConvectedDubins::Path * get_optPath(){return optPath_;};
  void get_optimalWaypoints( vd & x, vd & y, vd & hRad );
  void get_optimalWaypoints_XYEastNorth( vd & x, vd & y, vd & hRad );
  double get_optCost(){ return optPath_->get_T(); };


}; // class
} // namespace


#endif
