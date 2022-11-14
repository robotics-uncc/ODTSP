#include "ODTSP_Solver.h"
#include "ODTSP_Utils.h"


void ODTSP::Solver::writeSolnFile( std::string outputFilename ){
  FILE * fid = fopen(outputFilename.c_str(),"w");
  fprintf(fid,"solutionFileName='%s'\n",outputFilename.c_str());
  fprintf(fid,"algorithm=%d;\n",(int)alg_);
  
  if ( alg_ > 1 ){
    fprintf(fid,"tspSolver=%d;\n",(int)gtspSolver_);
  }
  else {
    fprintf(fid,"tspSolver=-1;\n");
  }

  fprintf(fid,"replanCostTotal=%3.3f;\n", path_.get_cost() );
  fprintf(fid,"gtspCostTotal=%3.3f;\n", gtspCost_ );
  fprintf(fid,"time=%3.3f;\n", timeToSolve_ );
  fprintf(fid,"listSize=%d;\n",list_.size());

  // path costs
  fprintf(fid,"pathCosts=[");
  vd pathCosts = path_.get_pathCosts();
  for ( int j = 0; j < pathCosts.size(); j++){
    fprintf(fid,"%3.3f,",pathCosts[j]);
  }
  fprintf(fid,"];\n");
  // orbit costs
  fprintf(fid,"orbitCosts=[");
  vd orbitCosts = path_.get_orbitCosts();
  for ( int j = 0; j < orbitCosts.size(); j++){
    fprintf(fid,"%3.3f,",orbitCosts[j]);
  }
  fprintf(fid,"];\n");

  for( int i = 0; i < list_.size(); i++ ){
    fprintf(fid,"list_elem{%d}=[",i+1);
    vi elems = list_.get_elems(i);
    for ( int j = 0; j < list_.get_numElems(i); j++){
      fprintf(fid,"%d,",elems[j]+1);
    }
    fprintf(fid,"];\n");
    // transfer x , y 
    fprintf(fid,"transferX{%d}=[",i+1);
    vd transferX = path_.get_xTransfer(i);
    for ( int j = 0; j < transferX.size(); j++){
      fprintf(fid,"%3.3f,",transferX[j]);
    }
    fprintf(fid,"];\n");    
    fprintf(fid,"transferY{%d}=[",i+1);
    vd transferY = path_.get_yTransfer(i);
    for ( int j = 0; j < transferY.size(); j++){
      fprintf(fid,"%3.3f,",transferY[j]);
    }
    fprintf(fid,"];\n");    
    // orbit x , y 
    fprintf(fid,"orbitX{%d}=[",i+1);
    vd orbitX = path_.get_xOrbit(i);
    for ( int j = 0; j < orbitX.size(); j++){
      fprintf(fid,"%3.3f,",orbitX[j]);
    }
    fprintf(fid,"];\n");
    fprintf(fid,"orbitY{%d}=[",i+1);
    vd orbitY = path_.get_yOrbit(i);
    for ( int j = 0; j < orbitY.size(); j++){
      fprintf(fid,"%3.3f,",orbitY[j]);
    }
    fprintf(fid,"];\n");
    fprintf(fid,"centerX(%d)=%3.3f;\n",i+1, list_.get_centerX(i));
    fprintf(fid,"centerY(%d)=%3.3f;\n",i+1, list_.get_centerY(i));
    fprintf(fid,"orbitRadius(%d)=%3.3f;\n",i+1, list_.get_orbitRadius(i));
    fprintf(fid,"orbitDir(%d)=%3.3f;\n",i+1, list_.get_orbitDir(i));
    fprintf(fid,"entryX(%d)=%3.3f;\n",i+1, list_.get_entryX(i) );
    fprintf(fid,"entryY(%d)=%3.3f;\n",i+1, list_.get_entryY(i) );
    fprintf(fid,"entryHrad(%d)=%3.3f;\n",i+1, list_.get_entryHrad(i) );
    fprintf(fid,"entryAngleRad(%d)=%3.3f;\n",i+1, list_.get_entryAngleRad(i));
    fprintf(fid,"exitX(%d)=%3.3f;\n",i+1, list_.get_exitX(i));
    fprintf(fid,"exitY(%d)=%3.3f;\n",i+1, list_.get_exitY(i));
    fprintf(fid,"exitHrad(%d)=%3.3f;\n",i+1, list_.get_exitHrad(i));
    fprintf(fid,"exitAngleRad(%d)=%3.3f;\n",i+1, list_.get_exitAngleRad(i));

  }
  fclose(fid);
}

bool ODTSP::Solver::processParams( ODTSP::RequiredParams params, 
                                   double orbitRadiusMax, double aspectAngleMaxRad ){
  prob_.set_targs( params.xTargs, params.yTargs );
  prob_.set_orbitProperties( params.orbitAngleRad, params.orbitRadius, orbitRadiusMax,    
                             aspectAngleMaxRad );
  filename_ = params.filename;
  precision_ = params.precision;
  return true;
}

bool ODTSP::Solver::processParams( ODTSP::RequiredParams params ){
  prob_.set_targs( params.xTargs, params.yTargs );
  prob_.set_orbitProperties( params.orbitAngleRad, params.orbitRadius );
  filename_ = params.filename;
  precision_ = params.precision;
  return true;
}



// SYM_ANGLES_ETSP
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::RequiredParams params ){
  processParams(params);
  if ( alg == ODTSP::Algorithm::SYM_ANGLES_ETSP ){
    alg_ = alg;
    prob_.set_motionProperties( params.minTurnRadius );    
    symAngleSolve();
  }
  else{
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
}

// SYM_ANGLES_ETSP_CURRENTS
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::RequiredParams params, 
                       double vehicleSpeed, double windMagnitude, double windDirRad ){
  processParams(params);
  if ( alg == ODTSP::Algorithm::SYM_ANGLES_ETSP_CURRENTS ){
    alg_ = alg;
    prob_.set_motionProperties( params.minTurnRadius, vehicleSpeed, ODTSP::motionModel::CONVECTED_DUBINS  );   
    prob_.set_windMagDir( windMagnitude, windDirRad ); 
    symAngleSolve();
  }
  else{
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
}

bool ODTSP::Solver::symAngleSolve(){
  ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
  ODTSP::ClusterMode clusterMode = ODTSP::ClusterMode::None;
  list_.initialize( &prob_, clusterMode);
  // solve etsp
  ODTSP::Dopt tspSolver( ODTSP::DoptSolvers::ETSP_NO_CLUSTERS, &prob_, &list_, precision_, filename_);
  std::vector<int> seq = tspSolver.get_atspSoln();    
  // use sequence to sort the targets 
  int N = prob_.get_numTargs();
  std::vector<double> xTargsSorted(N);
  std::vector<double> yTargsSorted(N);
  for (int i = 0; i < seq.size(); i++ ){
    xTargsSorted[i] = prob_.get_xTargs( seq[i] );
    yTargsSorted[i] = prob_.get_yTargs( seq[i] );
  }
  // analytical solution for entry angles    
  std::vector<double> entryAngles(N);
  std::vector<int> orbitDirections(N);
  symmetricEntryAnglesTour( xTargsSorted, yTargsSorted,
                            prob_.get_orbitAngleRad(), 
                            // output
                            entryAngles, orbitDirections);
  ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
  timeToSolve_ = (t1 - t0) / 1000000.0L;
  // set the list
  list_.set_orbitProperties( seq, entryAngles, orbitDirections);
  path_.set_problemList( &prob_ , &list_ );
  path_.set_saveWaypoints(true);
  path_.plan();
  path_.print();
  return true;
}

// NO_CLUSTER_GTSP
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver,
                       ODTSP::RequiredParams params, 
                       int numEntryAngles ){  
  if ( alg == ODTSP::Algorithm::NO_CLUSTER_GTSP ){
    alg_ = alg;
    gtspSolver_ = gtspSolver;
    prob_.set_motionProperties( params.minTurnRadius );
    processParams(params);
    ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
    ODTSP::ClusterMode clusterMode = ODTSP::ClusterMode::None;
    // set cluster mode
    list_.initialize( &prob_, clusterMode );
    // solve the gtsp
    ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      gtsp.solve( ODTSP::DoptSolvers::LKH_NO_CLUSTERS );
    }
    else {
      gtsp.solve( ODTSP::DoptSolvers::GLKH_NO_CLUSTERS );
    }
    ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
    timeToSolve_ = (t1 - t0) / 1000000.0L;
    // get solution outputs
    gtspCost_ = gtsp.computeCost();
    path_ = gtsp.get_path();      
  }
  else{
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
}

// NO_CLUSTER_GTSP_CURRENTS
// MAXIMAL_CLUSTER_TUNING
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
                       ODTSP::RequiredParams params, int numEntryAngles,
                       double val1, double val2, double val3 ){  
  // This function is overloaded as below with save number of arguments:
  //  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
  //          ODTSP::RequiredParams params, 
  //          int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad, 
  //          double tuningParam );
  //  Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
  //          ODTSP::RequiredParams params, int numEntryAngles,
  //          double vehicleSpeed, double windMagnitude, double windDirRad )

  if ( alg == ODTSP::Algorithm::NO_CLUSTER_GTSP ){
    double vehicleSpeed = val1;
    double windMagnitude = val2;
    double windDirRad = val3;
    processParams(params);
    alg_ = alg;
    gtspSolver_ = gtspSolver;
    // set cluster mode
    prob_.set_motionProperties( params.minTurnRadius, vehicleSpeed, ODTSP::motionModel::CONVECTED_DUBINS  );   
    prob_.set_windMagDir( windMagnitude, windDirRad ); 
    ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
    ODTSP::ClusterMode clusterMode = ODTSP::ClusterMode::None;
    list_.initialize( &prob_, clusterMode );
    // solve the gtsp
    ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      gtsp.solve( ODTSP::DoptSolvers::LKH_NO_CLUSTERS );
    }
    else {
      gtsp.solve( ODTSP::DoptSolvers::GLKH_NO_CLUSTERS );
    }
    ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
    timeToSolve_ = (t1 - t0) / 1000000.0L;
    // get solution outputs
    gtspCost_ = gtsp.computeCost();
    path_ = gtsp.get_path();    
  }
  else if ( alg == ODTSP::Algorithm::MAXIMAL_CLUSTER_TUNING_GTSP ){
    double orbitRadiusMax = val1;
    double aspectAngleMaxRad = val2;
    double tuningParam = val3;
    processParams( params, orbitRadiusMax, aspectAngleMaxRad );
    alg_ = alg;
    gtspSolver_ = gtspSolver;
    // set cluster mode
    prob_.set_motionProperties( params.minTurnRadius );   
    //prob_.set_windMagDir( windMagnitude, windDirRad ); 
    ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
    ODTSP::ClusterMode clusterMode = ODTSP::ClusterMode::MaximalClusterTuning;
    list_.set_tuningParam( tuningParam );
    list_.initialize( &prob_, clusterMode );
    // solve the gtsp
    ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      gtsp.solve( ODTSP::DoptSolvers::LKH_NO_CLUSTERS );
    }
    else {
      gtsp.solve( ODTSP::DoptSolvers::GLKH_NO_CLUSTERS );
    }
    ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
    timeToSolve_ = (t1 - t0) / 1000000.0L;
    // get solution outputs
    gtspCost_ = gtsp.computeCost();   
    path_ = gtsp.get_path();     
  }
  else {
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
}


// SIMPLE_CLUSTER_GTSP
// EXHAUSTIVE_CLUSTER_GTSP
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
                       ODTSP::RequiredParams params, 
                       int numEntryAngles, double orbitRadiusMax, double aspectAngleMaxRad ){  
  alg_ = alg;
  gtspSolver_ = gtspSolver;
  prob_.set_motionProperties( params.minTurnRadius );
  ODTSP::DoptSolvers doptSolverType;
  ODTSP::ClusterMode clusterMode;
  ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
  if ( alg == ODTSP::Algorithm::SIMPLE_CLUSTER_GTSP ){
    processParams( params, orbitRadiusMax, aspectAngleMaxRad );
    clusterMode = ODTSP::ClusterMode::SimpleCluster;    
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      doptSolverType = ODTSP::DoptSolvers::LKH_NO_CLUSTERS;
    }
    else {
      doptSolverType = ODTSP::DoptSolvers::GLKH_NO_CLUSTERS;
    }
  }
  else if ( alg == ODTSP::Algorithm::EXHAUSTIVE_CLUSTER_GTSP ){
    processParams( params, orbitRadiusMax, aspectAngleMaxRad );
    clusterMode = ODTSP::ClusterMode::Exhaustive;
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      doptSolverType = ODTSP::DoptSolvers::LKH_INTERSECTING_CLUSTERS;
    }
    else {
      doptSolverType = ODTSP::DoptSolvers::GLKH_INTERSECTING_CLUSTERS;
    }
  }
  else{
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
  // set cluster mode
  list_.initialize( &prob_, clusterMode );
  // solve the gtsp
  ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
  gtsp.solve( doptSolverType );
  ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
  timeToSolve_ = (t1 - t0) / 1000000.0L;
  // get solution outputs
  gtspCost_ = gtsp.computeCost();
  path_ = gtsp.get_path();  
  
}

// SIMPLE_CLUSTER_GTSP_CURRENTS
// EXHAUSTIVE_CLUSTER_GTSP_CURRENTS
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
                       ODTSP::RequiredParams params, int numEntryAngles,
                       double orbitRadiusMax, double aspectAngleMaxRad,
                       double vehicleSpeed, double windMagnitude, 
                       double windDirRad ){  
  alg_ = alg;
  gtspSolver_ = gtspSolver;
  processParams(params, orbitRadiusMax, aspectAngleMaxRad );
  // set cluster mode
  prob_.set_motionProperties( params.minTurnRadius, vehicleSpeed, ODTSP::motionModel::CONVECTED_DUBINS  );   
  prob_.set_windMagDir( windMagnitude, windDirRad ); 
  ODTSP::DoptSolvers doptSolverType;
  ODTSP::ClusterMode clusterMode;
  ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
  if ( alg == ODTSP::Algorithm::SIMPLE_CLUSTER_GTSP ){
    clusterMode = ODTSP::ClusterMode::SimpleCluster;    
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      doptSolverType = ODTSP::DoptSolvers::LKH_NO_CLUSTERS;
    }
    else {
      doptSolverType = ODTSP::DoptSolvers::GLKH_NO_CLUSTERS;
    }
  }
  else if ( alg == ODTSP::Algorithm::EXHAUSTIVE_CLUSTER_GTSP ){
    clusterMode = ODTSP::ClusterMode::Exhaustive;
    if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
      doptSolverType = ODTSP::DoptSolvers::LKH_INTERSECTING_CLUSTERS;
    }
    else {
      doptSolverType = ODTSP::DoptSolvers::GLKH_INTERSECTING_CLUSTERS;
    }
  }
  else{
    printf("ODTSP::Solver: Error, incorrect algorithm specified.\n");
  }
  // set cluster mode
  list_.initialize( &prob_, clusterMode );
  // solve the gtsp
  ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
  gtsp.solve( doptSolverType );
  ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
  timeToSolve_ = (t1 - t0) / 1000000.0L;
  // get solution outputs
  gtspCost_ = gtsp.computeCost();
  path_ = gtsp.get_path();    
}

// MAXIMAL_CLUSTER_TUNING_GTSP_CURRENTS
// MAXIMAL_CLUSTER_FEASIBLE_DUBINS_GTSP
ODTSP::Solver::Solver( ODTSP::Algorithm alg, ODTSP::GTSP_SOLVER gtspSolver, 
                       ODTSP::RequiredParams params, int numEntryAngles,
                       double orbitRadiusMax, double aspectAngleMaxRad,
                       double tuningParam, 
                       double vehicleSpeed, double windMagnitude, 
                       double windDirRad ){  
  alg_ = alg;
  gtspSolver_ = gtspSolver;
  double eps = windMagnitude / vehicleSpeed;
  double aspectAngleMaxRadMod = aspectAngleMaxRad - asin(eps);
  if ( aspectAngleMaxRadMod < 0 ){
    printf ("ODTSP::Solver::Solver : Currents too high problem not feasible.\n");
  }
  processParams(params, orbitRadiusMax, aspectAngleMaxRadMod );
  // set cluster mode
  if ( alg == ODTSP::Algorithm::MAXIMAL_CLUSTER_FEASIBLE_DUBINS_GTSP ){
    
    prob_.set_motionProperties( params.minTurnRadius, vehicleSpeed,                 
                                ODTSP::motionModel::FEASIBLE_DUBINS  );
  }
  else if ( alg == ODTSP::Algorithm::MAXIMAL_CLUSTER_TUNING_GTSP_CURRENTS ){
    prob_.set_motionProperties( params.minTurnRadius, vehicleSpeed,     
                                ODTSP::motionModel::CONVECTED_DUBINS  );   
  }
  prob_.set_windMagDir( windMagnitude, windDirRad ); 

  prob_.print();
  ODTSP::DoptSolvers doptSolverType;
  ODTSP::timestamp_t t0 = ODTSP::get_timestamp();
  ODTSP::ClusterMode clusterMode = ODTSP::ClusterMode::MaximalClusterTuning;
  list_.set_tuningParam( tuningParam );

  if ( gtspSolver == ODTSP::GTSP_SOLVER::LKH ){
    doptSolverType = ODTSP::DoptSolvers::LKH_NO_CLUSTERS;
  }
  else {
    doptSolverType = ODTSP::DoptSolvers::GLKH_NO_CLUSTERS;
  }
  // set cluster mode
  list_.initialize( &prob_, clusterMode );
  // solve the gtsp
  ODTSP::Dopt gtsp( &prob_, &list_, numEntryAngles, precision_, filename_ ); 
  gtsp.solve( doptSolverType );
  ODTSP::timestamp_t t1 = ODTSP::get_timestamp();
  timeToSolve_ = (t1 - t0) / 1000000.0L;
  // get solution outputs
  gtspCost_ = gtsp.computeCost();
  path_ = gtsp.get_path();    
}

