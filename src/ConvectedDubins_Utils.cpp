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

#include "ConvectedDubins_Utils.h"
#include "ConvectedDubins_Problem.h"
#include<cmath>
#include<stdexcept>


double ConvectedDubins::headingAngleFromInertialCourse( double v, double w, double hWRad, 
                                                        double chiRad){
  // see ODTSP paper 
  return fmodPos( chiRad - asin( w/v*sin(hWRad-chiRad) ), 2.0*M_PI );
}

double ConvectedDubins::speedOverGround( const double & v, 
                                         const double & w, 
                                         const double & hRad){
  // see ODTSP paper 
  double term = v*v + w*w + 2*v*w*cos(hRad);
  return sqrt(term*term);
}

double ConvectedDubins::courseAngleRad( const double & va, 
                                        const double & vw, 
                                        const double & hRad){
  // eq. 15 
  double num = va*sin(hRad);
  double denom = va*cos(hRad) + vw;
  return fmod( atan2(num, denom) , 2.0*M_PI);
}

void ConvectedDubins::convertInertialtoTrochoid( double xn, double ye,
                                                 double hwRad, 
                                                 double & xt, double & yt ){
    // from rotation matrix on p. 1737, above eq. 4
    xt =  cos(hwRad)*xn + sin(hwRad)*ye;
    yt = -sin(hwRad)*xn + cos(hwRad)*ye; 
}

double ConvectedDubins::trochoid_delx( double v, double omega, double wx, 
                                       double dir, double delt, double hInitRad ){
  return v/(dir*omega)*( sin(dir*omega*delt+hInitRad) - sin(hInitRad) )+wx*delt;
}

double ConvectedDubins::trochoid_delx( ConvectedDubins::Problem * prob, 
                                       double dir, double delt, double hInitRad ){
  return prob->v()/(dir*prob->omega())*( sin(dir*prob->omega()*delt+hInitRad) 
                                         - sin(hInitRad) ) + prob->wx()*delt;
}

double ConvectedDubins::trochoid_dely( double v, double omega, double wy,
                                       double dir, double delt, double hInitRad ){
  return v/(dir*omega)*(-cos(dir*omega*delt+hInitRad) + cos(hInitRad) )+wy*delt;
}

double ConvectedDubins::trochoid_dely( ConvectedDubins::Problem * prob, 
                                       double dir, double delt, double hInitRad ){
  return prob->v()/(dir*prob->omega())*( -cos(dir*prob->omega()*delt+hInitRad) + cos(hInitRad) )
                  + prob->wy()*delt;
}

double ConvectedDubins::trochoid_delh( double omega, 
                                       double dir, double delt ){
  return dir*(omega)*delt;
}

double ConvectedDubins::trochoid_delh( ConvectedDubins::Problem * prob, 
                                       double dir, double delt ){
  return dir*( prob->omega() )*delt;
}

double ConvectedDubins::straight_delx( double v, double wx, 
                                       double delt, double hInitRad ){
  return ( v*cos(hInitRad) + wx )*delt;
}

double ConvectedDubins::straight_delx( ConvectedDubins::Problem * prob, 
                                       double delt, double hInitRad ){
  return ( prob->v()*cos(hInitRad) + prob->wx() )*delt;
}

double ConvectedDubins::straight_dely( double v, double wy, 
                                       double delt, double hInitRad ){
  return ( v*sin(hInitRad) + wy )*delt;
}

double ConvectedDubins::straight_dely( ConvectedDubins::Problem * prob, 
                                       double delt, double hInitRad ){
  return ( prob->v()*sin(hInitRad) + prob->wy() )*delt;
}

ConvectedDubins::pathClass ConvectedDubins::getClassFromType( 
                                                  ConvectedDubins::pathType pt){
  switch (pt){
    case ConvectedDubins::pathType::LSL: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::LSR: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::RSL: return ConvectedDubins::pathClass::BSB; 
    case ConvectedDubins::pathType::RSR: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::LRL: return ConvectedDubins::pathClass::BBB;
    case ConvectedDubins::pathType::RLR: return ConvectedDubins::pathClass::BBB; 
  }

}

void ConvectedDubins::computeEndpoint( ConvectedDubins::Problem * prob, 
                                       ConvectedDubins::pathType pt,
                                       double dirA, double dirC,
                                       double tA, double tBeta, double T,
                                       double & xFinal, double & yFinal, 
                                       double & hFinalRad){

  // state after first trochoid straight
  double xBinit = prob->xInitial() + trochoid_delx( prob, dirA, tA, prob->hInitialRad() );
  double yBinit = prob->yInitial() + trochoid_dely( prob, dirA, tA, prob->hInitialRad() );
  double hBinit = prob->hInitialRad() + trochoid_delh( prob->omega(), dirA, tA ); 
  ConvectedDubins::pathClass pc = ConvectedDubins::getClassFromType(pt);


  double xBfinal, yBfinal, hBfinal;
  switch( pc ){
    // state after second segment (straight)
    case ConvectedDubins::pathClass::BSB: {
      xBfinal = xBinit + straight_delx( prob, tBeta - tA, hBinit );
      yBfinal = yBinit + straight_dely( prob, tBeta - tA, hBinit );
      hBfinal = hBinit;
      break;
    }
    // state after second segment (trochoid)
    case ConvectedDubins::pathClass::BBB: {
      xBfinal = xBinit + trochoid_delx( prob, -dirA, tBeta - tA , hBinit);
      yBfinal = yBinit + trochoid_dely( prob, -dirA, tBeta - tA , hBinit);
      hBfinal = hBinit + trochoid_delh( prob, -dirA, tBeta - tA ); 
      break;
    }
  }

  // final state 
  xFinal = xBfinal + trochoid_delx( prob, dirC, T-tBeta , hBfinal );
  yFinal = yBfinal + trochoid_dely( prob, dirC, T-tBeta , hBfinal );
  hFinalRad = hBfinal + trochoid_delh( prob, dirC, T-tBeta ); 
  hFinalRad = fmod(hFinalRad, 2.0*M_PI);
  if ( hFinalRad < 0 ){
    hFinalRad = hFinalRad + 2.0*M_PI;
  }
}

double ConvectedDubins::fmodPos( double val, double mod ){
  while ( val < 0 ){
    val = val + mod;
  }
  return fmod(val,mod);
}


