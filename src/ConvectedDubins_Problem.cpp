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

#include "ConvectedDubins_Problem.h"
#include "ConvectedDubins_Utils.h"
#include<assert.h>
#include<stdexcept>
#include<cmath>

void ConvectedDubins::Problem::initialize(){
  xInitial_ = 0.0;
  yInitial_ = 0.0;
  hInitialRad_ = 0.0;
  R_ = 1.0;
  v_ = 1.0;
  omega_ = R_/v_;
  w_ = 0.5;
  hwRad_ = M_PI/2.0;
  wx_ = w_*cos(hwRad_);
  wy_ = w_*sin(hwRad_);

  startPtInputFlag_ = false;
  endPtInputFlag_ = false;
  vehPropInputFlag_= false;
  envPropInputFlag_ = false;
}


void ConvectedDubins::Problem::set_xInitial( double xInitial ){
  xInitial_ = xInitial;
  startPtInputFlag_ = true;
}; 

void ConvectedDubins::Problem::set_yInitial( double yInitial ){
  yInitial_ = yInitial;
  startPtInputFlag_ = true;
};

void ConvectedDubins::Problem::set_hInitialRad( double hInitialRad ){
  #ifndef NDEBUG
  if ( (hInitialRad > 2.0*M_PI) || (hInitialRad < 0.0 ) ){
    throw std::runtime_error("ConvectedDubins::Problem::set_hInitialRad: Invalid heading.");
  }
  #endif
  hInitialRad_ = hInitialRad;
  startPtInputFlag_ = true;
};

void ConvectedDubins::Problem::set_xFinal( double xFinal ){
  xFinal_ = xFinal;
  endPtInputFlag_ = true;
}

void ConvectedDubins::Problem::set_yFinal( double yFinal ){
  yFinal_ = yFinal;
  endPtInputFlag_ = true;
}

void ConvectedDubins::Problem::set_hFinalRad( double hFinalRad ){
  #ifndef NDEBUG
  if ( (hFinalRad > 2.0*M_PI) || (hFinalRad < 0.0 ) ){
    throw std::runtime_error("ConvectedDubins::Problem::set_hFinalRad: Invalid heading.");
  }
  #endif
  hFinalRad_ = hFinalRad;
  endPtInputFlag_ = true;
}

void ConvectedDubins::Problem::printInputType( bool userInputFlag ){
  if ( userInputFlag ){
    printf(" [user input]\n");
  }
  else{
    printf(" [default]\n");
  }
}

ConvectedDubins::Problem::Problem(){
  initialize();
}

void ConvectedDubins::Problem::set_vehicleProperties( double minTurningRadius, 
                                                      double nominalSpeed ){
  #ifndef NDEBUG
  if ( minTurningRadius < 0 ){
    throw std::runtime_error("ConvectedDubins::Problem::sset_vehicleProperties: Invalid turn radius.");
  }
  if ( nominalSpeed <= 0){
    throw std::runtime_error("ConvectedDubins::Problem::set_nominalSpeed: Invalid speed.");
  }
  #endif 
	R_ = minTurningRadius;
  v_ = nominalSpeed;
  omega_ = v_/R_;
  vehPropInputFlag_ = true;
}

void ConvectedDubins::Problem::set_windMagDir( double windMagnitude, 
                                               double windDirRad ){
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
  w_ = windMagnitude;
  envPropInputFlag_ = true;
}

void ConvectedDubins::Problem::set_windMagXY( double windMagX,  
                                              double windMagY){
  w_ = std::sqrt(windMagX*windMagX + windMagY*windMagY);
  wx_ = windMagX;
  wy_ = windMagY;
  // TODO: hwRad_ = fmod( atan2(wx_,wy_) , 2.0*M_PI);   // 4-Jan-2019 changed:
  hwRad_ = fmod( atan2(wy_, wx_) , 2.0*M_PI);
  envPropInputFlag_ = true;
}

void ConvectedDubins::Problem::set_stateInitial(double xInitial, 
                                                double yInitial, 
                                                double hInitialRad ){
	set_xInitial(xInitial);
	set_yInitial(yInitial);
	set_hInitialRad(hInitialRad); 
}

void ConvectedDubins::Problem::set_stateInitial( const vd & stateInitial ){
	set_xInitial(stateInitial[0]);
	set_yInitial(stateInitial[1]);
	set_hInitialRad(stateInitial[2]);
}

void ConvectedDubins::Problem::set_stateInitial_XYEastNorth( double xInitial, 
                                                             double yInitial, 
                                                           double hInitialRad ){
	set_xInitial(yInitial);
	set_yInitial(xInitial);
	set_hInitialRad( fmodPos( M_PI/2.0-hInitialRad, 2.0*M_PI) ); 
}

void ConvectedDubins::Problem::set_stateInitial_XYEastNorthCourse( 
                                                         double xInitial, 
                                                         double yInitial, 
                                                         double chiInitialRad ){
  // check that wind has been input by user
  #ifndef NDEBUG
  if ( envPropInputFlag_ == false ){
    throw std::runtime_error("ConvectedDubins::Problem::set_stateInitial_XYEastNorthCourse: Set wind values before specifying course angle.");
  }
  #endif 

	set_xInitial(yInitial);
	set_yInitial(xInitial);
  // convert from ENU course angle to NED 
  double chiConverted = fmodPos( M_PI/2.0-chiInitialRad, 2.0*M_PI);
  // obtain heading in currents with this course 
  double hWithCurrents = headingAngleFromInertialCourse( v_, w_, hwRad_, 
                                                         chiConverted );
	set_hInitialRad( fmodPos( hWithCurrents , 2.0*M_PI) ); 
}

void ConvectedDubins::Problem::set_stateFinal( double xFinal, 
                                               double yFinal, 
                                               double hFinalRad ){
	set_xFinal(xFinal);
	set_yFinal(yFinal);
	set_hFinalRad(hFinalRad);
}

void ConvectedDubins::Problem::set_stateFinal( const vd & stateFinal ){
	set_xFinal(stateFinal[0]);
	set_yFinal(stateFinal[1]);
	set_hFinalRad(stateFinal[2]);
}

void ConvectedDubins::Problem::set_stateFinal_XYEastNorth( double xFinal, 
                                                           double yFinal, 
                                                           double hFinalRad ){
	set_xFinal(yFinal);
	set_yFinal(xFinal);
	set_hFinalRad( fmodPos( M_PI/2.0-hFinalRad, 2.0*M_PI) );
}

void ConvectedDubins::Problem::set_stateFinal_XYEastNorthCourse( double xFinal, 
                                                                 double yFinal, 
                                                           double chiFinalRad ){
	set_xFinal(yFinal);
	set_yFinal(xFinal);
  // see set_stateInitial_XYEastNorthCourse
  double chiConverted = fmodPos( M_PI/2.0-chiFinalRad, 2.0*M_PI);
  double hWithCurrents = headingAngleFromInertialCourse( v_, w_, hwRad_, chiConverted );
	set_hFinalRad( fmodPos( hWithCurrents , 2.0*M_PI) ); 
}

void ConvectedDubins::Problem::print(){
	printf("=========================================\n");
	printf("Convected Dubins Problem Statement \n");
	printf("=========================================\n");
  printf("(xInit, yInit, hInitDeg) : \t(%3.3f, %3.3f, %3.3f)",
            xInitial_, yInitial_, hInitialRad_*180/M_PI);
  printInputType(startPtInputFlag_);
  printf("(xFinal, yFinal, hFinalDeg) : \t(%3.3f, %3.3f, %3.3f)",xFinal_, 
            yFinal_, hFinalRad_*180/M_PI);
  printInputType(endPtInputFlag_);
  printf("Min. Turning Radius : \t\t%3.3f ", R_);
  printInputType(vehPropInputFlag_);
  printf("Vehicle Speed : \t\t%3.3f ", v_);
  printInputType(vehPropInputFlag_);
  printf("Max Turn Rate : \t\t%3.3f ", omega_);
  printInputType(vehPropInputFlag_);
  printf("Wind Speed : \t\t\t%3.3f ", w_);
  printInputType(envPropInputFlag_);
  printf("Wind Direction (deg.) : \t%3.3f ", hwRad_*180/M_PI);
  printInputType(envPropInputFlag_);
  printf("Wind Speed (wx,wy) : \t\t(%3.3f , %3.3f) ", wx_, wy_);
  printInputType(envPropInputFlag_);
}

bool ConvectedDubins::Problem::isDefined(){
  if ( !startPtInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: initial state not defined.\n");
    return false;
  }
  if ( !endPtInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: terminal state not defined.\n");
    return false;
  }
  if ( !vehPropInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: vehicle properties not defined.\n");
    return false;
  }
  if ( !envPropInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: environment properties not defined.\n");
    return false;
  }
  return true;
}

bool ConvectedDubins::Problem::isDefinedNoEndpoint(){
  if ( !startPtInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: initial state not defined.\n");
    return false;
  }
  if ( !vehPropInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: vehicle properties not defined.\n");
    return false;
  }
  if ( !envPropInputFlag_ ){
    printf("ConvectedDubins::Problem::isDefined: environment properties not defined.\n");
    return false;
  } 
  //printf(" return true \n");
  return true;
}


// get functions
vd ConvectedDubins::Problem::stateInitial(){
  vd stateInitial = {xInitial_, yInitial_, hInitialRad_};
  return stateInitial;
}

double ConvectedDubins::Problem::xFinal(){
  #ifndef NDEBUG
  if ( !endPtInputFlag_ ){
    throw std::runtime_error("ConvectedDubins::Problem::xFinal: Final state not defined.");
  }
  #endif
  return xFinal_;
};


double ConvectedDubins::Problem::yFinal(){
  #ifndef NDEBUG
  if ( !endPtInputFlag_ ){
    throw std::runtime_error("ConvectedDubins::Problem::yFinal: Final state not defined.");
  }
  #endif 
  return yFinal_;
};

double ConvectedDubins::Problem::hFinalRad(){
  #ifndef NDEBUG
  if ( !endPtInputFlag_ ){
    throw std::runtime_error("ConvectedDubins::Problem::hFinal: Final state not defined.");
  }
  #endif 
  return hFinalRad_;
};

vd ConvectedDubins::Problem::stateFinal(){
  #ifndef NDEBUG
  if ( !endPtInputFlag_ ){
    throw std::runtime_error("ConvectedDubins::Problem::stateFinal: Final state not defined.");
  }
  #endif 
  vd stateFinal = {xFinal_, yFinal_, hFinalRad_};
  return stateFinal;
}
