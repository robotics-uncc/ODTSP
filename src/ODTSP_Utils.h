// A. Wolek, 21-July-2016
// ODTSP_Utils.h

#ifndef ODTSP_UTILS_H
#define ODTSP_UTILS_H

// custom
#include "MathTools.h"
#include "ODTSP_Target.h"
#include "ODTSP_OrbitList.h"


#include <sys/time.h>


namespace ODTSP {

    typedef unsigned long long timestamp_t;


    // Usage:
    //    timestamp_t t0 = get_timestamp();
    //    // Process
    //    timestamp_t t1 = get_timestamp();
    //    double secs = (t1 - t0) / 1000000.0L;
    timestamp_t get_timestamp(void);

// * Note : LU = length unit

//// -----------------------------------------------------------------------------
//// transformEntryAngles: 
////    Converts given target locations and orbit entry angles to planar positions
////    and headings at the start and end of each orbit.
////
//// Inputs:    orbitRadius (LU) - radius of the orbit around each target
////            orbitAngle (rad) - angle of the arc subtended along each orbit
////            entryAngles(rad) - angle at which the orbit is started
////            orbitDirections (1/-1) - direction of orbit, 1 = CCW, -1 = CW
////            xTargs (LU)      - x-locations of the targets to be orbitted
////            yTargs (LU)      - y-locations of the targets to be orbitted
////  
//// Outputs:   entryX (LU)      - x-locations of the points where orbits start
////            entryY (LU)      - y-locations of the points where orbits start
////            entryH (rad)     - heading at which the orbits start
////            exitAngles (rad) - 
////            exitX (LU)       - x-locations of the points where orbits end
////            exitY (LU)       - y-locations of the points where orbits end
////            exitH (rad)      - heading at which the orbits end


//void transformEntryAngles(double orbitRadius, double orbitAngle,
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
//                                       std::vector<double> & exitH            );

//void transformEntryAngles(         const std::vector< ODTSP::Target> & targets, 
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
//                                   std::vector<double> & exitH            );

void transformEntryAngle(             const double & orbitRadius, 
                                      const double & entryAngle,
                                      const double & dir,
                                      const double & xTarg,
                                      const double & yTarg, 
                                      // outputs:
                                      double & orbitX,
                                      double & orbitY,
                                      double & orbitHrad);

//void transformEntryAngles( ODTSP::Problem * prob,
//                                  ODTSP::OrbitList * list );

//// -----------------------------------------------------------------------------
//// transformStatesToEntryAngles
////    Converts an unordered list of target locations, a target ordering,
////    and a list of points along the orbits corresponding to the ordering
////    into entry angles and an ordered target list
////
//// Inputs:    sequenceTargs     - integer list of targets to orbit 
////            xTargsList (LU)   - unordered list of x-locations of the targets
////            yTargsList (LU)   - unordered list of y-locations of the targets
////            fromX (LU)        - x-locations of the points where orbits exits
////            fromY (LU)        - y-locations of the points where orbits exits
////            fromH (rad)       - heading at which the orbits exits
////            orbitAngle (rad)  - orbit angle 
////  
//// Outputs:   entryAngles (rad)   - angle at which the orbit is started
////            orbitDirections (1/-1) - direction of orbit, 1 = CCW, -1 = CW
////            xTargsOrdered (LU)  - ordered list of x-locations of the targets
////            yTargsOrdered (LU)  - ordered list of y-locations of the targets
//void transformStatesToEntryAngles(std::vector<int>  & sequenceTargs,                                       
//                                  std::vector<double> & xTargsList,
//                                  std::vector<double> & yTargsList,
//                                  std::vector<double> & fromX,
//                                  std::vector<double> & fromY,
//                                  std::vector<double> & fromH,
//                                  double              & orbitAngle,
//                                  std::vector<double> & entryAngles,
//                                  std::vector<int>    & orbitDirections,
//                                  std::vector<double> & xTargsOrdered,
//                                  std::vector<double> & yTargsOrdered    );
//                                       

// -----------------------------------------------------------------------------
// symmetricEntryAngles:
//    Given a set of positions calculate the "symmetric" entry angles 
//
// Inputs:    orbitAngle (rad) - angle of the arc subtended along each orbit
//            xTargs (LU)      - x-locations of the targets to be orbitted
//            yTargs (LU)      - y-locations of the targets to be orbitted
//            xInit (LU)       - x-location of vehicle start position
//            yInit (LU)       - y-location of vehicle start position
//            xFinal (LU)      - x-location of vehicle end position
//            yFinal (LU)      - y-location of vehicle end position
//
// Outputs:   std::vector<double> - a vector of entryAngles(rad) angle at which 

// if start/end positions not supplied we assume this is a tour not a path
void symmetricEntryAnglesTour(       std::vector<double> xTargs,
                                     std::vector<double> yTargs,
                                     double orbitAngle,
                                     std::vector<double> & entryAngles,
                                     std::vector<int> & orbitDirections);

// the main function here computes the symmetric cangle for three points:
//
// Inputs:    Point A : (xA, yA)
//            Point B : (xB, yB)
//            Point C : (xC, yC)
//            orbitAngle (rad) - angle of the arc subtended along each orbit
//
// Output:    symmetrincEntryAngle (rad) angle at which the orbit is entered
//            so that the orbit is "symmetric" with respect to A-B-C
void symmetricEntryAngle(double xA, double yA, double xB, double yB,
                         double xC, double yC, double orbitAngle, 
                         double & entryAngle, int & orbitDirection );

double  symmetricEntryAngle(double xA, double yA, double xB, double yB,
                            double orbitAngle                             );

int dubinsPathSetToInteger(std::vector< std::string > dubinsPathSet);
int dubinsPathToInteger(std::string pathType);

vi nearestNeighborSeq( const vd & x, const vd & y);

double fmodPos( double val, double mod );


} // namespace

#endif

