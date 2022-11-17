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

#ifndef CONVECTED_DUBINS_UTILS_H
#define CONVECTED_DUBINS_UTILS_H

#include "ConvectedDubins_Path.h"
#include<vector>
#include<string>
#include<assert.h>

typedef std::vector<double> vd;
typedef std::vector<int> vi;

// The Utils are misc. utilities used by other functions throughout code 

namespace ConvectedDubins {

// Determines heading required for a give course angle:
//  v - flow-relative speed
//  w - disturbance (current/wind) speed
//  hwRad - direction of disturbance, measured CW from X axis pointing North 
//  chiRad - course angle (of vehicle's instantaenous motion)
double headingAngleFromInertialCourse( double v, double w, double hWRad, 
                                       double chiRad);

// Determines speed over ground at given heading angle 
double speedOverGround( const double & v, const double & w, 
                        const double & hRad);

// course angle in trochoidal reference frame (x-axis rotated downwind))
double courseAngleRad( const double & va, const double & vw, 
                       const double & hRad);

// convert (xn, ye) north-east coordinates to trochoidal (xt, yt) coordinates
// where hwRad is direction of wind 
//    delt is the amount of time spent along this extremal
void convertInertialtoTrochoid( double xn, double ye, double hwRad,
                                double & xt, double & yt );

// the following functions are used to compute the change in state in NED frame
// they are derived by integrating eq. 1 
double trochoid_delx ( double v, double omega, double wx, 
                       double dir, double delt, double hInitRad );

// same as above but uses (v, omega, wx) values contained in prob 
double trochoid_delx ( ConvectedDubins::Problem * prob, 
                       double dir, double delt, double hInitRad );

// see above 
double trochoid_dely ( double v, double omega, double wy, 
                       double dir, double delt, double hInitRad );
double trochoid_dely ( ConvectedDubins::Problem * prob, 
                       double dir, double delt, double hInitRad );

// similar but for heading change 
double trochoid_delh ( double omega, 
                       double dir, double delt );
double trochoid_delh ( ConvectedDubins::Problem * prob, 
                       double dir, double delt );

// similar but along straight segment, no heading change 
double straight_delx ( double v, double wx, 
                       double delt, double hInitRad );
double straight_delx ( ConvectedDubins::Problem * prob, 
                       double delt, double hInitRad );

double straight_dely ( double v, double wy, 
                       double delt, double hInitRad );
double straight_dely ( ConvectedDubins::Problem * prob, 
                       double delt, double hInitRad );

// returns a pathClass (BSB, BBB) given a path type (LSL, LSR, LRL, etc.)
ConvectedDubins::pathClass getClassFromType( ConvectedDubins::pathType pt );


// use the 'del' functions above to compute endpoint given parameters of a 
// convected dubins path. These parameters are:
//  pt : path type (LSL, LSR, etc.)
//  dirA : direction along first segment
//  dirC : direction along third segment 
//  tA : the end of first segment 
//  tBeta : the end of second segment
//  T : the end of third segment (total elapsed time))
//
// output is stored in arguments xFinal, yFinal, hFinal 
void computeEndpoint( ConvectedDubins::Problem * prob, 
                      ConvectedDubins::pathType pt,
                      double dirA, double dirC,
                      double tA, double tBeta, double T,
                      double & xFinal, double & yFinal, 
                      double & hFinalRad);

// like fmod (modulus), but ensures returned value is positive between [0, mod]
double fmodPos( double val, double mod );

} // namespace

#endif
