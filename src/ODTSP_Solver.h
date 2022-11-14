#include "ODTSP_Dopt.h"
#include "ODTSP_Problem.h"
#include "ODTSP_Path.h"
#include "ODTSP_Utils.h"


namespace ODTSP {

enum Algorithm {
  SYM_ANGLES_ETSP=0,
  SYM_ANGLES_ETSP_CURRENTS=1,
  NO_CLUSTER_GTSP=2,
  NO_CLUSTER_GTSP_CURRENTS=3,
  SIMPLE_CLUSTER_GTSP=4,
  SIMPLE_CLUSTER_GTSP_CURRENTS=5,
  EXHAUSTIVE_CLUSTER_GTSP=6,
  EXHAUSTIVE_CLUSTER_GTSP_CURRENTS=7,
  MAXIMAL_CLUSTER_TUNING_GTSP=8,
  MAXIMAL_CLUSTER_TUNING_GTSP_CURRENTS=9,
  MAXIMAL_CLUSTER_FEASIBLE_DUBINS_GTSP=10
};

enum GTSP_SOLVER {
  LKH=0,
  GLKH=1
};

struct RequiredParams {
  std::vector<double> xTargs;
  std::vector<double> yTargs;
  double orbitAngleRad;
  double orbitRadius;
  double minTurnRadius;
  std::string filename = "test";
  int precision = 1000;
};

class Solver{

public:


  // ************************ SOLVER OPTIONS **********************************
  // Note: You can all only call one of these functions

  // SYM_ANGLES_ETSP
  Solver( ODTSP::Algorithm alg, ODTSP::RequiredParams params );

  // SYM_ANGLES_ETSP_CURRENTS
  Solver( ODTSP::Algorithm alg, ODTSP::RequiredParams params, 
          double vehicleSpeed, double windMagnitude, double windDirRad );

  // NO_CLUSTER_GTSP
  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
          ODTSP::RequiredParams params, int numEntryAngles );

  // NO_CLUSTER_GTSP_CURRENTS
  // MAXIMAL_CLUSTER_TUNING
  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver,
          ODTSP::RequiredParams params, int numEntryAngles,
          double val1, double val2, double val3 ); 
  // This function is overloaded as below with save number of arguments:
  //  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
  //          ODTSP::RequiredParams params, 
  //          int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad, 
  //          double tuningParam );
  //  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
  //          ODTSP::RequiredParams params, int numEntryAngles,
  //          double vehicleSpeed, double windMagnitude, double windDirRad )

  // SIMPLE_CLUSTER_GTSP
  // EXHAUSTIVE_CLUSTER_GTSP
  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
          ODTSP::RequiredParams params, 
          int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad );

  // SIMPLE_CLUSTER_GTSP_CURRENTS
  // EXHAUSTIVE_CLUSTER_GTSP_CURRENTS
  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
          ODTSP::RequiredParams params, 
          int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad,
          double vehicleSpeed, double windMagnitude, double windDirRad );

  // MAXIMAL_CLUSTER_TUNING_CURRENTS
  // MAXIMAL_CLUSTER_FEASIBLE_DUBINS_GTSP
  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
          ODTSP::RequiredParams params, 
          int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad,
          double tuningParam, 
          double vehicleSpeed, double windMagnitude, double windDirRad );


  ~Solver(){};


  // ************************ OUTPUT OPTIONS **********************************

  // this output will give waypoints 
  ODTSP::Path* getPath(){return &path_;};
  ODTSP::OrbitList* getOrbitList(){return &list_;};
  void writeSolnFile( std::string outputFilename );


private:
  
  ODTSP::Algorithm alg_;
  ODTSP::GTSP_SOLVER gtspSolver_;
  ODTSP::Problem prob_;
  ODTSP::Path path_;
  ODTSP::OrbitList list_;
 
  double gtspCost_ = -1.0;

  // sub-solvers
  bool symAngleSolve(); // SYM_ANGLES_ETSP and SYM_ANGLES_ETSP_CURRENTS


  // no clustering 
  bool processParams( ODTSP::RequiredParams params );
  // with clustering 
  bool processParams( ODTSP::RequiredParams params, 
                      double orbitRadiusMax, double aspectAngleMaxRad ); 

  // for gtsp
  std::string filename_;
  int precision_;  

  double timeToSolve_;
  


};

} //namespace

