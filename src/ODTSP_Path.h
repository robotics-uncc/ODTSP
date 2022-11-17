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

// A. Wolek, 29-Aug-2016
// ODTSP_Path.h

#ifndef ODTSP_PATH_H
#define ODTSP_PATH_H

#include "ODTSP_Problem.h"
#include "ODTSP_Target.h"
#include "ODTSP_OrbitList.h"

#include "ConvectedDubins_Problem.h"
#include "ConvectedDubins_Solver.h"
#include "ConvectedDubins_Path.h"

#include "RobustDubins_Problem.h"
#include "RobustDubins_Solver.h"
#include "RobustDubins_Path.h"
#include "MathTools.h"

namespace ODTSP {

class Path {

private:

  // inputs
  ODTSP::Problem* prob_; // the ODTSP problem
  ODTSP::OrbitList* list_; // gives the set of targets, sequence, entry angles, 
                           // and orbit directions
  
  void initialize();

  // optional:
  bool saveWaypoints_;   // flag indicating if the waypoints should be computed
  int numWpts_; // the number of waypoints for each orbit (circular arc)
  bool saveDubinsPaths_; // flag indicating if m_DubinsPaths should
                          // be accessible after plan() is executed

  // given a ODTSP::Problem and entry angles there is a corresponding set of 
  // Dubins Problems defined by the following vector. 
  // The main purpose of plan() is to construct and solve these sub-problems. 
  std::vector<RobustDubins::Problem> problemSet_; 
  std::vector<ConvectedDubins::Problem> problemSetConvected_; 

  // intermediate planning functions
  void planTour(); // used when init/final state not specified for ODTSProblem

  // plans the path along a particular orbit  
  void planOrbit(const double & xTarg, const double & yTarg, 
                 const double & entryAngle, const double & exitAngle, 
                 const double & orbitRadius, const int & dir);

  // plans a Dubins Path connecting two orbits
  void planDubinsPath(RobustDubins::Problem & dubinsProb);
  void planCDubinsPath(ConvectedDubins::Problem & cDubinsProb);

  // outputs
  std::vector<double> pathCosts_; // cost of intermediate Dubins paths
  std::vector<double> orbitCosts_; // cost of intermediate Dubins paths
  double              totalCost_; // summation of all Dubins path costs

  // optional outputs
  // (x,y,h) states along the path
  std::vector<double> xOverallPath_, yOverallPath_, hOverallPath_;
 
 std::vector< std::vector<double> > transferX_;
 std::vector< std::vector<double> > transferY_;
 std::vector< std::vector<double> > orbitX_;
 std::vector< std::vector<double> > orbitY_;

  
  bool pathPlanned_;

  // corresponding vector of Dubins path objects
  std::vector< RobustDubins::Path > dubinsPaths_;
  std::vector< std::string  >       dubinsPathTypes_;

  std::vector< ConvectedDubins::Path* > cDubinsPaths_;
public: 

  // constructor, destructor
  Path(ODTSP::Problem* prob, ODTSP::OrbitList* list);
  Path();
  ~Path(){};
  // required if not using 1-st form of constructor
  void set_problemList(ODTSP::Problem* prob, ODTSP::OrbitList* list);

  // optional set functions
  void set_numWpts(int numWpts){numWpts_ = numWpts;};
  void set_saveWaypoints(bool flag){saveWaypoints_ = flag;};
  void set_saveDubinsPaths(bool flag){saveDubinsPaths_ = flag;};

  // main process
  void plan();
  void print();
  void writeCommands(std::string filename);
  void writePlotCommands(std::string filename);

  // get inputs
  ODTSP::Problem*     get_ODTSProblem(){return prob_;};
  ODTSP::OrbitList*   get_orbitList(){return list_;};
  int                 get_numWpts(){return numWpts_;};

  // this returns cost of orbit transfers (not including orbit)
  std::vector<double> get_pathCosts(){return pathCosts_;};
  std::vector<double> get_orbitCosts(){return orbitCosts_;};

  // this returns sum(pathCosts_)
  double              get_cost(){return totalCost_;};

  // get optional outputs
  std::vector<RobustDubins::Problem> get_problemSet(){return problemSet_;};
  std::vector<double> get_xOverallPath(){return xOverallPath_;};
  std::vector<double> get_yOverallPath(){return yOverallPath_;};
  std::vector<double> get_hOverallPath(){return hOverallPath_;};

  std::vector<double> get_xTransfer(int ind){return transferX_[ind];};
  std::vector<double> get_yTransfer(int ind){return transferY_[ind];};

  std::vector<double> get_xOrbit(int ind){return orbitX_[ind];};
  std::vector<double> get_yOrbit(int ind){return orbitY_[ind];};

  std::vector< RobustDubins::Path > get_DubinsPaths(){return dubinsPaths_;};
  std::vector< std::string  > get_DubinsPathTypes(){
    return dubinsPathTypes_;
  };
  int get_DubinsPathTypesID();

};

} // namespace

#endif

