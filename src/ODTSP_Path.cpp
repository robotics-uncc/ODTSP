#include "ODTSP_Path.h"
#include "ODTSP_Utils.h"
#include "ODTSP_Problem.h"

#include <stdio.h>

ODTSP::Path::Path(ODTSP::Problem* prob, ODTSP::OrbitList* list){
  initialize();
  // intialize based on user input
  prob_ = prob;
  list_ = list;
}

ODTSP::Path::Path(){
  initialize();
}

void ODTSP::Path::set_problemList(ODTSP::Problem* prob, ODTSP::OrbitList* list){
  prob_ = prob;
  list_ = list;
}

void ODTSP::Path::initialize(){
  // initialize default values
  numWpts_ = 100;
  saveWaypoints_ = false;
  saveDubinsPaths_ = false;
  problemSet_.clear();
  pathCosts_.clear();
  orbitCosts_.clear();
  transferX_.clear();
  transferY_.clear();
  orbitX_.clear();
  orbitY_.clear();
  totalCost_ = -1;
  xOverallPath_.clear(); 
  yOverallPath_.clear();
  hOverallPath_.clear();
  dubinsPaths_.clear(); 
  dubinsPathTypes_.clear();
}

void ODTSP::Path::plan(){
  planTour();
  totalCost_ = MathTools::sum(pathCosts_) + MathTools::sum(orbitCosts_); // calculate total cost
  pathPlanned_ = true;
}

void ODTSP::Path::planTour(){
  printf("*************** PLAN TOUR ***************\n");
  // plan inter-target transfers
  for (int i = 0; i < list_->size()-1; i++){
    if ( prob_->motionModel() == DUBINS ){
      RobustDubins::Problem dubinsProb; // define the problem
      dubinsProb.set_minTurningRadius(prob_->R());
      // update problem statement
      dubinsProb.set_stateInitial( list_->get_exitX(i), 
                                   list_->get_exitY(i), 
                                   list_->get_exitHrad(i) );
      dubinsProb.set_stateFinal( list_->get_entryX(i+1), 
                                 list_->get_entryY(i+1), 
                                 list_->get_entryHrad(i+1) );
      planDubinsPath(dubinsProb);
    }
    else if ( prob_->motionModel() == FEASIBLE_DUBINS ){
      RobustDubins::Problem dubinsProb; // define the problem
      double eps = prob_->w() / prob_->v();
      double feasibleRadius = prob_->R() * (1.0 + eps) * (1.0 + eps);
      dubinsProb.set_minTurningRadius( feasibleRadius );
      // update problem statement
      dubinsProb.set_stateInitial( list_->get_exitX(i), 
                                   list_->get_exitY(i), 
                                   list_->get_exitHrad(i) );
      dubinsProb.set_stateFinal( list_->get_entryX(i+1), 
                                 list_->get_entryY(i+1), 
                                 list_->get_entryHrad(i+1) );
      planDubinsPath(dubinsProb);
    }
    else if ( prob_->motionModel() == CONVECTED_DUBINS ){
      ConvectedDubins::Problem cDubinsProb;
      // Warning!!! list_->get_exitHrad and list_->get_entryHrad are actually
      // course angles on the circle. 
      cDubinsProb.set_vehicleProperties( prob_->R() , prob_->v() );
      cDubinsProb.set_windMagDir( prob_->w() , prob_->hwRadNED() );
      cDubinsProb.set_stateInitial_XYEastNorthCourse( list_->get_exitX(i), 
                                                      list_->get_exitY(i), 
                                                      list_->get_exitHrad(i) );
      cDubinsProb.set_stateFinal_XYEastNorthCourse( list_->get_entryX(i+1), 
                                                    list_->get_entryY(i+1), 
                                                    list_->get_entryHrad(i+1) );
      planCDubinsPath(cDubinsProb);
    }
    planOrbit( list_->get_centerX(i+1), list_->get_centerY(i+1), 
               list_->get_entryAngleRad(i+1), list_->get_exitAngleRad(i+1), 
               list_->get_orbitRadius(i+1), list_->get_orbitDir(i+1) );

  } 
  if ( prob_->motionModel() == DUBINS ){
    RobustDubins::Problem dubinsProb; // define the problem
    dubinsProb.set_minTurningRadius(prob_->R());
    // plan path from final target to intial target entry point
    dubinsProb.set_stateInitial( list_->get_exitX(list_->size()-1), 
                                 list_->get_exitY(list_->size()-1), 
                                 list_->get_exitHrad(list_->size()-1) );
    dubinsProb.set_stateFinal( list_->get_entryX(0), 
                               list_->get_entryY(0), 
                               list_->get_entryHrad(0) );
    planDubinsPath(dubinsProb);
  }
  else if ( prob_->motionModel() == FEASIBLE_DUBINS ){
    RobustDubins::Problem dubinsProb; // define the problem
    double eps = prob_->w() / prob_->v();
    double feasibleRadius = prob_->R() * (1.0 + eps) * (1.0 + eps);
    dubinsProb.set_minTurningRadius( feasibleRadius );
    // plan path from final target to intial target entry point
    dubinsProb.set_stateInitial( list_->get_exitX(list_->size()-1), 
                                 list_->get_exitY(list_->size()-1), 
                                 list_->get_exitHrad(list_->size()-1) );
    dubinsProb.set_stateFinal( list_->get_entryX(0), 
                               list_->get_entryY(0), 
                               list_->get_entryHrad(0) );
    planDubinsPath(dubinsProb);
  }
  else if ( prob_->motionModel() == CONVECTED_DUBINS ){
      ConvectedDubins::Problem cDubinsProb;
      cDubinsProb.set_vehicleProperties( prob_->R() , prob_->v() );
      cDubinsProb.set_windMagDir( prob_->w() , prob_->hwRadNED() );
      cDubinsProb.set_stateInitial_XYEastNorthCourse( list_->get_exitX(list_->size()-1), 
                                                      list_->get_exitY(list_->size()-1), 
                                                      list_->get_exitHrad(list_->size()-1) );
      cDubinsProb.set_stateFinal_XYEastNorthCourse( list_->get_entryX(0), 
                                                    list_->get_entryY(0), 
                                                    list_->get_entryHrad(0) );
      planCDubinsPath(cDubinsProb);
  }
  planOrbit( list_->get_centerX(0), list_->get_centerY(0), 
             list_->get_entryAngleRad(0), list_->get_exitAngleRad(0), 
             list_->get_orbitRadius(0), list_->get_orbitDir(0) );  
}

void ODTSP::Path::planCDubinsPath(ConvectedDubins::Problem & cDubinsProb){
  ConvectedDubins::Solver cds( &cDubinsProb );
  pathCosts_.push_back( cds.get_optCost() ); // store path cost
  cds.printSolns();
  // save the waypoints, if requested
  if ( saveWaypoints_ ){
    // get waypoints/states of the current path
    std::vector<double> xCurrentPath, yCurrentPath, hCurrentPath;
    cds.get_optimalWaypoints_XYEastNorth(xCurrentPath, yCurrentPath, hCurrentPath); 
    // append to overall path
    transferX_.push_back(xCurrentPath);
    transferY_.push_back(yCurrentPath);

    MathTools::append(xOverallPath_,xCurrentPath);
    MathTools::append(yOverallPath_,yCurrentPath);
    MathTools::append(hOverallPath_,hCurrentPath);
  }
  // save a copy of the entire dubins path object, if desired
  if ( saveDubinsPaths_ ){
    cDubinsPaths_.push_back( cds.get_optPath() );
  }
}

void ODTSP::Path::planDubinsPath(RobustDubins::Problem & dubinsProb){
  RobustDubins::Solver rds; // create the solver
  rds.set_problemStatement(dubinsProb); 
  rds.solve();
  if ( prob_->motionModel() == DUBINS ){
  pathCosts_.push_back( rds.get_optCost() ); // store path cost
  }
  else if ( prob_->motionModel() == FEASIBLE_DUBINS ){
  // input for optPathTimeInCurrents requires current in ENU not NED 
  pathCosts_.push_back( rds.get_optPathTimeInCurrents( prob_->v(), prob_->w(), prob_->hwRad() ) ); // store path cost
  }
  // save the waypoints, if requested
  if ( saveWaypoints_ ){
    // get waypoints/states of the current path
    std::vector<double> xCurrentPath, yCurrentPath, hCurrentPath;
    rds.get_optimalWaypoints(xCurrentPath, yCurrentPath, hCurrentPath); 
    // append to overall path
    transferX_.push_back(xCurrentPath);
    transferY_.push_back(yCurrentPath);

    MathTools::append(xOverallPath_,xCurrentPath);
    MathTools::append(yOverallPath_,yCurrentPath);
    MathTools::append(hOverallPath_,hCurrentPath);
  }
  // save a copy of the entire dubins path object, if desired
  if ( saveDubinsPaths_ ){
    dubinsPaths_.push_back( rds.get_optimalPath() );
    dubinsPathTypes_.push_back( (dubinsPaths_.at(
                                    dubinsPaths_.size() - 1)).get_pathType() );
  }
}

void ODTSP::Path::planOrbit(const double & xTarg, const double & yTarg, 
                            const double & entryAngle,
                            const double & exitAngle, 
                            const double & orbitRadius, const int & dir){
  if (saveWaypoints_){
    // get the parameters defining an arc
    std::vector<double> arcParams = {xTarg, yTarg, orbitRadius,
                                     entryAngle, exitAngle};

    // this function creates a set of angles using polarspace spanning from 
    // entryAngle to exitAngle using the direction dir 
    //   - it assumes the interval requested is less than one revolution 
    //   - 'dir' if +1 is CCW if -1 is CW 
    //   - it adds +/- 2pi to the final angle to achieve continuity 

    std::vector< std::vector<double> > arcWpts = MathTools::generateCircularArc(
                                                                arcParams, 
                                                                numWpts_, dir);
    MathTools::append(xOverallPath_,arcWpts.at(0));
    MathTools::append(yOverallPath_,arcWpts.at(1));
    // TODO: Append to hOverallPath_

    orbitX_.push_back( arcWpts.at(0) );
    orbitY_.push_back( arcWpts.at(1) );

  }
  // TODO: add cost of orbits here, check if dubins of c dubins 

  double cost = -1;
  
  // define the Dubins problem and solve
  if ( prob_->motionModel() == DUBINS){
    // then orbit cost is length 
    double orbitLength = orbitRadius*prob_->get_orbitAngleRad();
    cost = ( orbitLength );
  }
  else if ( prob_->motionModel() == CONVECTED_DUBINS || prob_->motionModel() == FEASIBLE_DUBINS ){

    // get the "from" point (i-th cluster, k-th node) + orbitAngle
    // (the x,y coordiantes are NOT scaled up)

    double dir2 = (double) dir;
    double toX, toY, toHrad;
    ODTSP::transformEntryAngle(orbitRadius, 
                               entryAngle, dir2, 
                               xTarg, 
                               yTarg, 
                               // output
                               toX, toY, toHrad);

    // not that the above transformEntryAngle/orbitlist function uses ENU frame
    // convection of X - East , angles +CCW from East 
    // dir +1 implies CCW 

    // but the below computation derived in paper assumes NED frame 
    // convention of X - North, angles +CW from North
    // dir +1 implies CW 

    // convert to XY = NE coordinate system 
    double v = prob_->v();
    double eps = prob_->w()/v; 
    double chiNED = fmod( M_PI/2.0 - toHrad , 2.0*M_PI );
    double varphiNED = prob_->hwRadNED(); // psi is measured from X axis north CW
    double psiNED;

    
    // convert dir
    double dirMod = -dir;    
    int N = 100;
    double delth = prob_->get_orbitAngleRad()/( double(N) );
    double timeSum = 0;
    
      
    double vg, delt;
    for ( int i = 0; i < N; i ++ ){

      psiNED = fmod( chiNED - asin( eps*sin(varphiNED - chiNED) ), 2.0*M_PI );
      vg = v*sqrt(1.0 + eps*( eps + 2.0*cos(varphiNED - psiNED) ) ) ;   
      delt = orbitRadius/vg*delth;
      // update 
      timeSum = timeSum + delt;
      chiNED = chiNED + delth*dirMod;
    }
    cost = timeSum;
  }
  if ( cost < 0 ){
    throw std::runtime_error("ODTSP::Path::planOrbit < 0\n"); 
  }
  orbitCosts_.push_back(cost);



}


// -----------------------------------------------------------------------------
void ODTSP::Path::print(){
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Dubins Tour Circle RID Path:" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << " Entry and Exit Angles (deg.) : " << std::endl;
  for (int i = 0; i < list_->size(); i++){
  std::cout << "(entry" << i << "," << "exit" << i << ") = (" 
            << list_->get_entryAngleRad(i)*MathTools::RAD2DEG
            << ", " << list_->get_exitAngleRad(i)*MathTools::RAD2DEG << ")"<< std::endl;
  }
  if ( saveDubinsPaths_){ 
    if ( prob_->motionModel() == DUBINS ){ 
      std::cout << " Path Init States \n (x \t y \t h, deg): " << std::endl;
      for (int i = 0; i < list_->size(); i++){
        printf("%6.6f \t %6.6f \t %6.6f \n",(dubinsPaths_.at(i)).get_xInitial()
               ,(dubinsPaths_.at(i)).get_yInitial()
               ,(dubinsPaths_.at(i)).get_hInitial()*MathTools::RAD2DEG); 
      }
      std::cout << " Path Final States \n (x \t y \t h, deg): " << std::endl;
      for (int i = 0; i < list_->size(); i++){
        printf("%6.6f \t %6.6f \t %6.6f \n",(dubinsPaths_.at(i)).get_xFinal()
               ,(dubinsPaths_.at(i)).get_yFinal()
               ,(dubinsPaths_.at(i)).get_hFinal()*MathTools::RAD2DEG); 
      }
      for (int i = 0; i < dubinsPaths_.size(); i++){
      std::cout <<">>> Dubins Path No. " << i << std::endl;
      std::cout <<"\t Type: " << (dubinsPaths_.at(i)).get_pathType()<< std::endl;
      std::cout <<"\t Cost: " << (dubinsPaths_.at(i)).get_cost()<< std::endl;
      printf(" Initial (x,y,hdeg) = (%3.3f, %3.3f, %3.3f) \n",
             (dubinsPaths_[i]).get_xInitial(), (dubinsPaths_[i]).get_yInitial(),
             (dubinsPaths_[i]).get_hInitial() * 180/M_PI  );
      printf(" Final (x,y,hdeg) = (%3.3f, %3.3f, %3.3f) \n",
             (dubinsPaths_[i]).get_xFinal(), (dubinsPaths_[i]).get_yFinal(),
             (dubinsPaths_[i]).get_hFinal() * 180/M_PI  );
      std::cout <<"\t a: "<<(dubinsPaths_.at(i)).get_aParamUnsigned()<<std::endl;
      std::cout <<"\t b: "<<(dubinsPaths_.at(i)).get_bParamUnsigned()<<std::endl;
      std::cout <<"\t c: "<<(dubinsPaths_.at(i)).get_cParamUnsigned()<<std::endl;
      }  
    }
  }
  std::cout << "Total Tour Cost : " 
            << totalCost_ << std::endl;
  std::cout << "---------------------------------------------" << std::endl;  
}


//// -----------------------------------------------------------------------------
//void ODTSP::Path::writeCommands(std::string filename){
//  // open file for writing
//  FILE * pFile; 
//  pFile = fopen (filename.c_str(),"w");
//  for (int i = 0; i < dubinsPaths_.size(); i++){
//    std::cout << " i : " << i << std::endl;
//    // write orbit command
//    std::string dir = "CCW";
//    if (orbitDirections_.at(i) == -1){
//      dir = "CW";
//    } 
//    // ORBIT, ORBIT_CENTER_X, ORBIT_CENTER_Y, ENTRYANGLE, DIR, ANGLE(ABS)
//    fprintf(pFile, "ORBIT_%d,%6.6f,%6.6f,%6.6f,%s,%6.6f\n",i, 
//                                               (targets_[i]).xc_, 
//                                               (targets_[i]).yc_,
//                                               entryAngles_.at(i), dir.c_str(),
//                                               prob_->get_orbitAngle() );
//    // write dubins command
//    std::vector<double> segX, segY, segH;
//    std::vector<double> orbitCenterX, orbitCenterY, entryAngles, orbitDir;
//    (dubinsPaths_.at(i)).get_segPts(segX, segY, segH);
//    (dubinsPaths_.at(i)).get_orbitParams(orbitCenterX, orbitCenterY, 
//                                                         entryAngles, orbitDir);
//    std::string dpType = (dubinsPaths_.at(i)).get_pathType();

//    int turnCount = 1;
//    if ( (dubinsPaths_.at(i)).get_aParamUnsigned() > 0 ){
//      dir = "CCW";
//      if (orbitDir.at(0) == -1){
//        dir = "CW";
//      } 
//      fprintf(pFile, "DPTURN_%d_SEG%d,%6.6f,%6.6f,%6.6f,%s,%6.6f\n", i, turnCount, 
//                                     orbitCenterX.at(0),  orbitCenterY.at(0), 
//                                     entryAngles.at(0), dir.c_str(),
//                                     (dubinsPaths_.at(i)).get_aParamUnsigned());
//      turnCount = turnCount + 1;
//    }
//    if ( (dubinsPaths_.at(i)).get_bParamUnsigned() > 0 ){
//      if ( ((dubinsPaths_.at(i)).get_pathType()).compare("BBB")==0 ){
//        dir = "CCW";
//        if (orbitDir.at(1) == -1){
//          dir = "CW";
//        } 
//        fprintf(pFile, "DPTURN_%d_SEG%d,%6.6f,%6.6f,%6.6f,%s,%6.6f\n", i, turnCount, 
//                                       orbitCenterX.at(1), orbitCenterY.at(1), 
//                                       entryAngles.at(1), dir.c_str(),
//                                     (dubinsPaths_.at(i)).get_bParamUnsigned());
//        turnCount = turnCount + 1;
//      }
//      else {
//        fprintf(pFile, "DPSTRAIGHT_%d,%6.6f,%6.6f,%6.6f,%6.6f\n", i, segX.at(1), 
//                                                       segY.at(1), segX.at(2), 
//                                                       segY.at(2) );
//      }
//    }
//    if ( (dubinsPaths_.at(i)).get_cParamUnsigned() > 0 ){
//      int j = 1;
//      if ( ((dubinsPaths_.at(i)).get_pathType()).compare("BBB")==0 ){
//        j = 2;
//      }
//      dir = "CCW";
//      if (orbitDir.at(j) == -1){
//        dir = "CW";
//      } 
//      std::cout << " j " << j << std::endl;
////      (dubinsPaths_.at(i)).print();
//      fprintf(pFile, "DPTURN_%d_SEG%d,%6.6f,%6.6f,%6.6,%s,%6.6f\n", i, 
//                                    turnCount, orbitCenterX.at(j),  
//                                    orbitCenterY.at(j), 
//                                    entryAngles.at(j), dir.c_str(),
//                                    (dubinsPaths_.at(i)).get_cParamUnsigned());
//    }
//  }
//}
// -----------------------------------------------------------------------------
int ODTSP::Path::get_DubinsPathTypesID(){
  return ODTSP::dubinsPathSetToInteger(dubinsPathTypes_);
}
// -----------------------------------------------------------------------------
void ODTSP::Path::writePlotCommands(std::string fileName){
  if ( pathPlanned_ && saveWaypoints_ ){
    printf(">>>> writing plot commands... \n");
    prob_->print();
    std::vector<double> xData = prob_->get_xTargs();
    std::vector<double> yData = prob_->get_yTargs();
    for (int i = 0 ; i < xData.size(); i++ ){
      printf("(x,y)=(%3.3f,%3.3f)\n",xData[i],yData[i]);
    }
    MathTools::plot1DCurve(fileName, 1, prob_->get_xTargs(), 
                           prob_->get_yTargs(), "xt", "yt", "ro");
    // draw circles for each orbit
    for (int i = 0; i < list_->size(); i++){
      std::vector<double> circleParams = { list_->get_centerX(i), 
                                           list_->get_centerY(i), 
                                           list_->get_orbitRadius(i) };
      std::vector< std::vector<double> > circle = MathTools::generateCircle(
                                                                     circleParams, 
                                                                     100 );
      std::vector<double> circle_xpts = circle.at(0);
      std::vector<double> circle_ypts = circle.at(1);
      MathTools::plot1DCurve(fileName, 1, circle_xpts, circle_ypts, "xc", "yc", 
                                                                          "r--");      
    }

    MathTools::plot1DCurve(fileName, 1, xOverallPath_, yOverallPath_,"x", "y", 
                           "k-");
    MathTools::axisProperty(fileName, "equal");
  }
  else {
    if ( !pathPlanned_ ){
      throw std::runtime_error("Path::writePlotCommands: Path has not been planned yet: try running plan() first");
    }
    if ( !saveWaypoints_ ){
      throw std::runtime_error("Path::writePlotCommands: Path waypoints not saved: try set_saveWaypoints(true) before executing plan()");
    }
  }
}

