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

#include "ODTSP_Problem.h"
#include "ODTSP_Utils.h"

#include<MathTools.h>
#include<iostream>
#include<assert.h>

ODTSP::Problem::Problem(){
  numTargs_ = 0;
  xTargs_.clear();
  yTargs_.clear();
  orbitAngleRad_= 0;
  minTurnRadius_ = 0;
  orbitRadius_ = 0;
  orbitRadiusMax_ = 0;
  aspectAngleMaxRad_ = 0;
  dmax_ = 0;
  mm_ = UNDEFINED_MOTION_MODEL;
  v_ = 0;
  w_ = -1;
}

bool ODTSP::Problem::isDefined(){
  if ( numTargs_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: numTargs\n");
    return false;
  }
  if ( orbitAngleRad_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: orbitAngleRad\n");
    return false;
  }
  if ( orbitRadius_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: orbitRadius\n");
    return false;
  }
  if ( orbitRadiusMax_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: orbitRadiusMax\n");
    return true; // valid if clustering not used 
  }
  if ( aspectAngleMaxRad_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: aspectAngleMaxRad\n");
    return true; // valid if clustering not used 
  }
  if ( mm_ == UNDEFINED_MOTION_MODEL ){
    printf("ODTSP::Problem: Warning, problem ill-defined: UNDEFINED_MOTION_MODEL\n");
    return false;
  }
  if ( minTurnRadius_ == 0 ){
    printf("ODTSP::Problem: Warning, problem ill-defined: minTurnRadius\n");
    return false;
  }
  if ( mm_ == CONVECTED_DUBINS){
    if ( v_ == 0 ){
      printf("ODTSP::Problem: Warning, problem ill-defined: CONVECTED_DUBINS (v = 0)\n");
      return false;
    }
    if ( w_ == -1 ){
      printf("ODTSP::Problem: Warning, problem ill-defined: CONVECTED_DUBINS (w = 0)\n");
      return false;
    }
  }
  return true;
}

// set functions 

bool ODTSP::Problem::set_targs( const vd & xTargs, const vd & yTargs){
  #ifndef NDEBUG
  if ( xTargs.size() != yTargs.size() ){
    throw std::runtime_error("Error, Problem::set_targs: inputs are inconsistent size\n");
  }
  #endif
  xTargs_ = xTargs;
  yTargs_ = yTargs;
  numTargs_ = xTargs.size();
  return true;
}

bool ODTSP::Problem::set_orbitProperties( double orbitAngleRad, double orbitRadius,
                                          double orbitRadiusMax, double aspectAngleMaxRad ){
  #ifndef NDEBUG
  if ( orbitAngleRad <= 0 ){
    throw std::runtime_error("Error, Problem::set_orbitProperties: orbitAngleDeg <= 0\n");
  }
  if ( orbitRadiusMax < orbitRadius ){
    throw std::runtime_error("Error, Problem::set_orbitProperties: orbitRadiusMax < orbitRadius\n");
  } 
  if ( aspectAngleMaxRad > M_PI/2.0  || aspectAngleMaxRad < 0 ){  
    printf(" aspectAngleMaxRad = %3.3f rad (%3.3f deg) \n", aspectAngleMaxRad, aspectAngleMaxRad * 180.0 / M_PI );
    throw std::runtime_error("Error, Problem::set_orbitProperties: aspectAngleMaxDeg > pi/2 or < 0\n");
  }  
  #endif
  orbitAngleRad_ = orbitAngleRad;
  orbitRadius_ = orbitRadius;
  aspectAngleMaxRad_ = aspectAngleMaxRad;
  orbitRadiusMax_ = orbitRadiusMax;
  // computes the critical (dmax) seperation distance for joining targets into a
  // cluster. If the target centers are displaced by a distace d > dmax then 
  // they can never be joined into a cluster and can be treated as isolated pts.
  double dcrit1 = orbitRadiusMax_*( sin(aspectAngleMaxRad_) / 
                   (1+sin(aspectAngleMaxRad_)) ); // aspect constraint
  double dcrit2 = (orbitRadiusMax_ - orbitRadius_) / 2; // range constraint  
  dmax_ = std::min(dcrit1, dcrit2);
  return true;
}


bool ODTSP::Problem::set_orbitProperties( double orbitAngleRad, double orbitRadius  ){
  #ifndef NDEBUG
  if ( orbitAngleRad <= 0 ){
    throw std::runtime_error("Error, Problem::set_orbitProperties: orbitAngleDeg <= 0\n");
  }
  #endif
  orbitAngleRad_ = orbitAngleRad;
  orbitRadius_ = orbitRadius;
  return true;
}

bool ODTSP::Problem::set_motionProperties( double minTurnRadius ){
  #ifndef NDEBUG
  if ( minTurnRadius <= 0 ){
    throw std::runtime_error("Error, Problem::set_motionProperties: minTurnRadius <= 0\n");
  }
  #endif
  minTurnRadius_ = minTurnRadius;
  mm_ = ODTSP::DUBINS;
  return true;
}

bool ODTSP::Problem::set_motionProperties( double minTurnRadius, double vehicleSpeed,  ODTSP::motionModel mm ){
  #ifndef NDEBUG
  if ( minTurnRadius <= 0 ){
    throw std::runtime_error("Error, Problem::set_motionProperties: minTurnRadius <= 0\n");
  }
  if ( vehicleSpeed <= 0 ){
    throw std::runtime_error("Error, Problem::set_motionProperties: vehicleSpeed <= 0\n");
  }
  #endif
  minTurnRadius_ = minTurnRadius;
  v_ = vehicleSpeed;
  mm_ = mm;
  return true;
}


void ODTSP::Problem::set_windMagDir( double windMagnitude, double windDirRad ){
  #ifndef NDEBUG
  if ( windMagnitude < 0 ){
    throw std::runtime_error("ConvectedDubins::Problem::set_windMagDir: Invalid magnitude.");
  }
  if ( (windDirRad > 2.0*M_PI) || (windDirRad < 0.0 ) ){
    throw std::runtime_error("ConvectedDubins::Problem::set_windMagDir: Invalid direction.");
  }
  #endif 
  wx_ = windMagnitude*cos(windDirRad);
  wy_ = windMagnitude*sin(windDirRad);
  hwRad_ = windDirRad;
  hwRadNED_ = fmodPos( M_PI/2.0-hwRad_, 2.0*M_PI);
  w_ = windMagnitude;
  return;
}

void ODTSP::Problem::set_windMagXY( double windMagX, double windMagY){
  w_ = std::sqrt(windMagX*windMagX + windMagY*windMagY);
  wx_ = windMagX;
  wy_ = windMagY;
  hwRad_ = fmod( atan2(wx_,wy_) , 2.0*M_PI); 
  hwRadNED_ = fmodPos( M_PI/2.0-hwRad_, 2.0*M_PI); 
  return;
}


double  ODTSP::Problem::get_orbitRadiusMax(){
  #ifndef NDEBUG
  if ( orbitRadiusMax_ == 0 ){
    throw std::runtime_error("ODTSP::Problem::get_orbitRadiusMax: max orbit radius not set!");
  }
  #endif
  return orbitRadiusMax_;
}

double  ODTSP::Problem::get_aspectAngleMaxRad(){
  #ifndef NDEBUG
  if ( aspectAngleMaxRad_ == 0 ){
    throw std::runtime_error("ODTSP::Problem::get_aspectAngleMaxRad: aspect angle max not set!");
  }
  #endif
  return aspectAngleMaxRad_;
}

double  ODTSP::Problem::get_dmax(){
  #ifndef NDEBUG
  if ( aspectAngleMaxRad_ == 0 ){
    throw std::runtime_error("ODTSP::Problem::get_dmax: max orbit / aspect angle not set!");
  }
  #endif
  return dmax_;
}

void ODTSP::Problem::print(){
  printf("ODTSP Problem Statement:\n");
  printf("---------------------------------------------\n");
  printf("Number of Targets: %i\n", numTargs_);
  for (int i = 0; i < numTargs_; i++){
    printf("(x%d,y%d) = (%3.3f,%3.3f)\n",i,i,xTargs_[i],yTargs_[i]);
  }
  printf("Orbit Angle (deg): %3.3f\n", orbitAngleRad_*180.0/M_PI );
  printf("Orbit Radius (min.): %3.3f\n", orbitRadius_);
  if ( orbitRadiusMax_ > 0 ){
    // clustering is being considered with these properties
    printf("Orbit Radius (max.): %3.3f\n",orbitRadiusMax_);
    printf("Aspect Angle (max. deg): %3.3f\n",aspectAngleMaxRad_*180.0/M_PI);
    printf("Critical Distance: %3.3f\n", dmax_);
  }
   switch ( mm_ ){
      case(CONVECTED_DUBINS): {
        printf("Motion Model : Convected Dubins \n"); 
        printf("Vehicle Min Turn Radius: %3.3f\n", minTurnRadius_);
        printf("Vehicle Speed: %3.3f\n", v_);
        printf("Wind Speed and Direction (deg): (%3.3f , %3.3f) \n", w_, hwRad_*180.0/M_PI );
      } break;
      case(DUBINS): {
        printf("Motion Model: Dubins \n"); 
        printf("Vehicle Min Turn Radius: %3.3f\n", minTurnRadius_);      
      } break;
   } 
  return;
}
