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

// A. Wolek, 31-Aug-2016
// ODTSP_Utils.cpp

#include "ODTSP_Utils.h"
#include<algorithm>    // std::sort
#include<math.h>       /* pow */


ODTSP::timestamp_t ODTSP::get_timestamp(void){
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }


//// -----------------------------------------------------------------------------
//// transformEntryAngles
//void ODTSP::transformEntryAngles(      double orbitRadius, 
//                                       double orbitAngle,
//                                       std::vector<double> entryAngles,
//                                       std::vector<int>    orbitDirections,
//                                       std::vector<double> xTargs,
//                                       std::vector<double> yTargs,
//                                       std::vector<double> & entryX,
//                                       std::vector<double> & entryY,
//                                       std::vector<double> & entryH,
//                                       std::vector<double> & exitAngles,
//                                       std::vector<double> & exitX,
//                                       std::vector<double> & exitY,
//                                       std::vector<double> & exitH            ){
//  int numOrbits = entryAngles.size(); // number of targets
//  for (int i = 0; i < numOrbits; i++){  
//    double dir = orbitDirections.at(i);
//    exitAngles.at(i) = ( MathTools::mod(entryAngles.at(i) + dir*orbitAngle, 
//                                                                     2.0*M_PI));
//    entryX.at(i) = ( xTargs.at(i) + orbitRadius*cos(entryAngles.at(i)) );
//    entryY.at(i) = ( yTargs.at(i) + orbitRadius*sin(entryAngles.at(i)) );
//    entryH.at(i) = ( MathTools::mod(entryAngles.at(i)+dir*M_PI/2.0,2.0*M_PI) );
//    exitX.at(i) = ( xTargs.at(i) + orbitRadius*cos(exitAngles.at(i)) );
//    exitY.at(i) = ( yTargs.at(i) + orbitRadius*sin(exitAngles.at(i)) );
//    exitH.at(i) = ( MathTools::mod(exitAngles.at(i) + dir*M_PI/2.0, 2.0*M_PI));
//  }
//}


//// -----------------------------------------------------------------------------
//// transformEntryAngles
//void ODTSP::transformEntryAngles( ODTSP::Problem * prob,
//                                  ODTSP::OrbitList * list ){
//  // inputs  
//  std::vector<double> centerX = list->get_centerX();
//  std::vector<double> centerY = list->get_centerY();
//  std::vector<double> entryAngles = list->get_entryAngles();
//  std::vector<int> orbitDir = list->get_orbitDir();
//  std::vector<double> orbitRadius = list->get_orbitRadius();
//  int numOrbits = list->get_numOrbits();
//  double orbitAngle = prob->get_orbitAngle();
//  // outputs
//  std::vector<double> entryX(numOrbits);
//  std::vector<double> entryY(numOrbits);
//  std::vector<double> entryHrad(numOrbits);

//  std::vector<double> exitAngles(numOrbits);
//  std::vector<double> exitX(numOrbits);
//  std::vector<double> exitY(numOrbits);
//  std::vector<double> exitHrad(numOrbits);
//  // compute
//  for (int i = 0; i < numOrbits; i++){  
//    double dir = (double)orbitDir[i];
//    exitAngles[i] = ( MathTools::mod(entryAngles[i] + dir*orbitAngle,2.0*M_PI));
//    entryX[i] = ( centerX[i] + orbitRadius[i]*cos( entryAngles[i] ) );
//    entryY[i] = ( centerY[i] + orbitRadius[i]*sin( entryAngles[i] ) );
//    entryHrad[i] = ( MathTools::mod(entryAngles.at(i)+dir*M_PI/2.0,2.0*M_PI) );
//    exitX[i] = ( centerX[i] + orbitRadius[i]*cos( exitAngles[i] ) );
//    exitY[i] = ( centerY[i] + orbitRadius[i]*sin( exitAngles[i] ) );
//    exitHrad[i] = ( MathTools::mod( exitAngles[i] + dir*M_PI/2.0, 2.0*M_PI));
//  }
//  // copy back to orbit list
//  list->set_entryX(entryX);
//  list->set_entryY(entryY);
//  list->set_entryHrad(entryHrad);
//  list->set_exitX(exitX);
//  list->set_exitY(exitY);
//  list->set_exitHrad(exitHrad);
//  list->set_exitAngles(exitAngles);
//  list->set_statesPopulated(true);
//}

//// -----------------------------------------------------------------------------
//// transformEntryAngles
//void ODTSP::transformEntryAngles(  const std::vector< ODTSP::Target> & targets, 
//                                   const double & orbitAngle,
//                                   const std::vector<double> & entryAngles,
//                                   const std::vector<int>    & orbitDirections,
//                                   // outputs
//                                   std::vector<double> & entryX,
//                                   std::vector<double> & entryY,
//                                   std::vector<double> & entryH,
//                                   std::vector<double> & exitAngles,
//                                   std::vector<double> & exitX,
//                                   std::vector<double> & exitY,
//                                   std::vector<double> & exitH            ){
//  int numOrbits = entryAngles.size(); 
//  for (int i = 0; i < numOrbits; i++){  
//    double dir = orbitDirections.at(i);
//    exitAngles.at(i) = ( MathTools::mod(entryAngles.at(i) + dir*orbitAngle, 
//                                                                     2.0*M_PI));
//    entryX.at(i) = ((targets[i]).xc_ + (targets[i]).rc_*cos(entryAngles.at(i)));
//    entryY.at(i) = ((targets[i]).yc_ + (targets[i]).rc_*sin(entryAngles.at(i)));
//    entryH.at(i) = (MathTools::mod(entryAngles.at(i)+dir*M_PI/2.0,2.0*M_PI) );
//    exitX.at(i) = ((targets[i]).xc_ + (targets[i]).rc_*cos(exitAngles.at(i)) );
//    exitY.at(i) = ((targets[i]).yc_ + (targets[i]).rc_*sin(exitAngles.at(i)) );
//    exitH.at(i) = (MathTools::mod(exitAngles.at(i) + dir*M_PI/2.0, 2.0*M_PI));
//  }
//}


// -----------------------------------------------------------------------------
// transformEntryAngles
void ODTSP::transformEntryAngle(      const double & orbitRadius, 
                                      const double & entryAngle,
                                      const double & dir,
                                      const double & xTarg,
                                      const double & yTarg, 
                                      // outputs:
                                      double & orbitX,
                                      double & orbitY,
                                      double & orbitHrad){
    // entryAngle is assumed to be in radians, measured CCW from X axis East
    orbitX = ( xTarg + orbitRadius*cos(entryAngle) );
    orbitY = ( yTarg + orbitRadius*sin(entryAngle) );
    // orbitHrad gives heading, measured CCW from X axis East
    // dir is (+1) CCW, 
    // (e.g., dir = 1 and entryAngle = 0, implies heading is pi/2 (north) )
    orbitHrad = ( MathTools::mod(entryAngle+dir*M_PI/2.0, 2.0*M_PI) );
}



//// -----------------------------------------------------------------------------
//// transformStatesToEntryAngles
//void ODTSP::transformStatesToEntryAngles(std::vector<int>  & sequenceTargs,
//                                       std::vector<double> & xTargsList,
//                                       std::vector<double> & yTargsList,
//                                       std::vector<double> & fromX,
//                                       std::vector<double> & fromY,
//                                       std::vector<double> & fromH,
//                                       double              & orbitAngle,
//                                       std::vector<double> & entryAngles,
//                                       std::vector<int>    & orbitDirections,
//                                       std::vector<double> & xTargsOrdered,
//                                       std::vector<double> & yTargsOrdered){
//  // iterate through the given target ordering        
//  for (int i = 0; i < sequenceTargs.size(); i++){
//    int curTarg = sequenceTargs.at(i); // current target in the ordering
//    xTargsOrdered.at(i) = (xTargsList.at(curTarg)); // add to the ordered list
//    yTargsOrdered.at(i) = (yTargsList.at(curTarg));
//    // compute the entry angle measured from the current target center to the 
//    // (x,y) position along theorbit
//    double dely = fromY.at(i)-yTargsOrdered.at(i);
//    double delx = fromX.at(i)-xTargsOrdered.at(i);
//    double exitAngle = atan2( dely , delx );
//    // consider two possible directions of travel
//    // the sign of the entry angle determiens if the heading corresponds to a 
//    // clockwise (CW) or counter-clockwise (CCW) orbit
//    double dir = 1; // default (CCW)
//    double tol = 0.01; // account for round-off error in entryH 
//    // case 1: motion is CCW, the corresponding heading is:
//    double exitH_CCW = exitAngle + dir*M_PI/2.0;
//    // check if this is consistent with the entryH supplied
//    if ( MathTools::polarDistance(exitH_CCW, fromH.at(i)) > tol ){
//      // case 2: motion is CW
//      dir = -1; // if not consistent then the heading must be CW
//    }
//    entryAngles.at(i) = (MathTools::mod(exitAngle - dir*orbitAngle,2.0*M_PI));
//    orbitDirections.at(i) = (dir);
//  }
//}


// -----------------------------------------------------------------------------
// symmetricEntryAnglesTour
void ODTSP::symmetricEntryAnglesTour(std::vector<double> xTargs,
                                     std::vector<double> yTargs,
                                     double orbitAngle,
                                     std::vector<double> & entryAngles,
                                     std::vector<int> & orbitDirections){  
  int numOrbits = xTargs.size(); // number of targets 
  assert( numOrbits == yTargs.size() );
  assert( numOrbits == entryAngles.size() ); // must already be correct size
  assert( numOrbits == orbitDirections.size() );
  if ( numOrbits < 2 ){
    std::cout << "Error: insufficient number of orbits specified" << std::endl;
  }  
  else if ( numOrbits == 2 ){
    entryAngles.at(0)=( ODTSP::symmetricEntryAngle(xTargs.at(0),yTargs.at(0), 
                                                   xTargs.at(1), yTargs.at(1), 
                                                   orbitAngle ) );
    orbitDirections.at(0) = -1;
    entryAngles.at(1)=( ODTSP::symmetricEntryAngle(xTargs.at(1),yTargs.at(1), 
                                                   xTargs.at(0), yTargs.at(0), 
                                                   orbitAngle ) );
    orbitDirections.at(1) = -1;
  }
  else if ( numOrbits > 2 ){
    printf(" numOrbits = %d \n", numOrbits );
    double curEntryAngle;
    int curOrbitDirection;
    // append first target to end of vector to close the loop 
    xTargs.insert(xTargs.begin(), xTargs.at(xTargs.size()-1) );  
    yTargs.insert(yTargs.begin(), yTargs.at(yTargs.size()-1) ); 
    xTargs.insert(xTargs.end(), xTargs.at(1) );  
    yTargs.insert(yTargs.end(), yTargs.at(1) ); 
    for (int i = 1; i <=numOrbits; i++){ // iterate for each target triplet
      ODTSP::symmetricEntryAngle(xTargs.at(i-1), yTargs.at(i-1), xTargs.at(i), 
                                 yTargs.at(i), xTargs.at(i+1), yTargs.at(i+1), 
                                 orbitAngle, curEntryAngle, curOrbitDirection);
      entryAngles.at(i-1) = curEntryAngle;
      orbitDirections.at(i-1) = curOrbitDirection;
    }
  }
  return;
}
// -----------------------------------------------------------------------------
// symmetricEntryAngles
void ODTSP::symmetricEntryAngle(double xA, double yA, double xB, double yB,
                                double xC, double yC, double orbitAngle,
                                double & entryAngle, int & orbitDirection){  
    // suppose every triplet of target is labeled A-B-C. 
    // let thetaBA be the angle from B to A, measured from the x-axis
    // let thetaBC be the angle from B to C, measured from the x-axis
    double thetaBA = atan2(yA - yB, xA - xB);
    double thetaBC = atan2(yC - yB, xC - xB);
    thetaBC = MathTools::mod(thetaBC, 2.0*M_PI);
    thetaBA = MathTools::mod(thetaBA, 2.0*M_PI);
    // let alpha be the angle formed by ABC
    double alpha = std::max(thetaBA, thetaBC) - std::min(thetaBA, thetaBC);
    // calculate the entry angle so that the orbit is symmetric about the 
    // the ray that bisects ABC
    double bisectorAngle = std::min(thetaBA, thetaBC) + alpha/2.0;
    bisectorAngle = MathTools::mod(bisectorAngle, 2.0*M_PI);
    // the bisector should be on "inside" the interval of angles given by ABC
    if ( MathTools::polarDistance(bisectorAngle, thetaBA) > 
         MathTools::polarDistance(bisectorAngle + M_PI, thetaBA) ){
      bisectorAngle = MathTools::mod(bisectorAngle + M_PI, 2.0*M_PI);
    }
    double gamma = MathTools::mod( thetaBC + ( 2.0*M_PI - thetaBA ) , 2.0*M_PI);
    if ( gamma <= M_PI ){
      orbitDirection = -1;
      entryAngle = MathTools::mod( bisectorAngle+M_PI+orbitAngle/2 , 2.0*M_PI);
    }
    else {
      orbitDirection = 1;
      entryAngle = MathTools::mod( bisectorAngle-M_PI-orbitAngle/2 , 2.0*M_PI);
    }
    return;
}
// -----------------------------------------------------------------------------
// symmetricEntryAngles
double ODTSP::symmetricEntryAngle(double xA, double yA, double xB, double yB,
                                  double orbitAngle){  
    // let thetaBA be the angle from B to A, measured from the x-axis
    double thetaBA = atan2(yA - yB, xA - xB);
    thetaBA = MathTools::mod(thetaBA, 2.0*M_PI);
    // calculate entry angle
    double entryAngle = thetaBA - orbitAngle/2;
    return MathTools::mod(entryAngle, 2.0*M_PI);
}

// -----------------------------------------------------------------------------
// dubinsPathPairsToInteger
int ODTSP::dubinsPathSetToInteger(std::vector< std::string > dubinsPathSet){
  // initialize variables
  std::vector<int> pathVals;
  // get the path values 
  for (int i = 0; i < dubinsPathSet.size(); i++){
    pathVals.push_back(ODTSP::dubinsPathToInteger(dubinsPathSet.at(i))); 
  }
  // sort from ascending order
  std::sort(pathVals.begin(), pathVals.end());
  // assign a unique ID
  int uniqueID = 0;
  for (int i = 0; i < pathVals.size(); i++){
    uniqueID = pathVals.at(i)*pow(10.0, i) + uniqueID;
  }
  return uniqueID; 
}

int ODTSP::dubinsPathToInteger(std::string pathType){
  if ( pathType.compare("LSL")==0 ){return 1;}
  else if ( pathType.compare("LSR")==0 ){return 2;}
  else if ( pathType.compare("RSL")==0 ){return 3;}
  else if ( pathType.compare("RSR")==0 ){return 4;}
  else if ( pathType.compare("LRL")==0 ){return 5;}
  else if ( pathType.compare("RLR")==0 ){return 6;}
}


vi ODTSP::nearestNeighborSeq( const vd & x, const vd & y){
  assert ( x.size() == y.size() );  // error checking 
  int n = x.size();  // size of the x, y list 
  vi seq( n ); // intialize the sequence to be returned 
  vi targs( n ); // create a vector of targets (not including the first targ 0) 
  for (int i = 1; i < n; i++ ){
    targs[i] = i;
  }
  // determine the nearest neighbor sequence 
  seq[0] = 0;
  int k = 0;
  while ( targs.size() > 0 ){
    double minD = std::numeric_limits<double>::max();
    int minDInd = -1;
    for ( int j = 0; j < targs.size(); j++ ){
      double d = (x[seq[k]] - x[targs[j]])*(x[seq[k]] - x[targs[j]]) 
                 + (y[seq[k]] - y[targs[j]])*(y[seq[k]] - y[targs[j]]);
      if ( d < minD ){
        minD = d;
        minDInd = j;
      }
    }  
    seq[k] = targs[minDInd];
    targs.erase(targs.begin()+minDInd);
    k++;
  }
  return seq;
}

double ODTSP::fmodPos( double val, double mod ){
  while ( val < 0 ){
    val = val + mod;
  }
  return fmod(val,mod);
}
