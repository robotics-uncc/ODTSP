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

#include "ODTSP_Target.h"
#include<iostream>
#include<algorithm>

void ODTSP::Target::initialize(const std::vector<int> & elems, const double & centerX, 
                      const double & centerY, const double & orbitRadius,
                      const double & orbitAngleRad, 
                      const ODTSP::TargetType & type){
  // check if constructor arguments are valid
  #ifndef NDEBUG
  if ( elems.size() == 0 ){
    throw std::runtime_error("Error, Target Construction: element list of size zero!");
  }
  if ( orbitRadius < 0 ){
    throw std::runtime_error("Error, Target Construction: orbit radius negative!");
  }
  if ( orbitAngleRad < 0 ){
    throw std::runtime_error("Error, Target Construction: orbit angel negative!");
  }
  #endif
  // make a private copy
  elems_ = elems;
  centerX_ = centerX;
  centerY_ = centerY;
  orbitRadius_ = orbitRadius;
  orbitAngleRad_ = orbitAngleRad;
  type_ = type;
  // initialize default values
  tolAngleRad_ = 0.01*2.0*M_PI;
  tolLength_ = 0.01*orbitRadius_;
  targetDefined_ = false; 
}

bool ODTSP::Target::set_tolerance(const double & tolAngleRad, 
                                  const double & tolLength){
  #ifndef NDEBUG
  // do not allow negative ang. tolerance or tolerance greater than 5 deg.
  if ( tolAngleRad_ < 0 || tolAngleRad_ > 5.0*M_PI/180.0 ){
    char msg[150]; 
    sprintf(msg, "Error, Target::setTolerance: tolAngleRad invalid, got : %i", tolAngleRad_);
    throw std::runtime_error(msg);
  }  
  // do not allow negative length tol. or tol. greater than 1/8 the orbitRadius
  if ( tolLength < 0 || tolLength > 1/8*orbitRadius_){
    printf("Error, Target::setTolerance: tolLength invalid\n");
    char msg[150]; 
    sprintf(msg, "Error, Target::setTolerance: tolLength invalid, got : %i", tolLength);
    throw std::runtime_error(msg);
  }
  #endif
  // update default values
  tolAngleRad_ = tolAngleRad;
  tolLength_ = tolLength;
  return true;
}

bool ODTSP::Target::set_entryAngleOrbitDir(const double & entryAngleRad, 
                                           const int & orbitDir){
  // check if heading is a valid value
  #ifndef NDEBUG
  if ( !validAngle(entryAngleRad) ){
    char msg[150]; 
    sprintf(msg, "Error, Target::set_entryAngleOrbitDir: must be in radians on the interval [0,2*pi], but got : %f", entryAngleRad);
    throw std::runtime_error(msg);
  }  
  if ( orbitDir != 1 && orbitDir != -1 ){
    char msg[150]; 
    sprintf(msg, "Error, Target::set_entryAngleOrbitDir: orbitDir invalid, must be either 1 or -1, but got : %i", orbitDir);
    throw std::runtime_error(msg);
  }
  #endif
  // make a private copy
  orbitDir_ = orbitDir;
  entryAngleRad_ = entryAngleRad;
  // define target states
  exitAngleRad_ = MathTools::mod( entryAngleRad_ + (double)orbitDir*orbitAngleRad_, 2.0*M_PI ); 
  entryX_ = centerX_ + orbitRadius_*cos( entryAngleRad_ );
  entryY_ = centerY_ + orbitRadius_*sin( entryAngleRad_ );
  entryHrad_ = ( MathTools::mod( entryAngleRad_ 
                            + (double)orbitDir*M_PI/2.0, 2.0*M_PI) ); 
  exitX_ = centerX_ + orbitRadius_*cos( exitAngleRad_ );
  exitY_ = centerY_ + orbitRadius_*sin( exitAngleRad_ );
  exitHrad_ = ( MathTools::mod( exitAngleRad_ 
                            + (double)orbitDir*M_PI/2.0, 2.0*M_PI) );   
  targetDefined_ = true;
  return true; 
}

bool ODTSP::Target::set_entryState(const double & entryX, 
                                   const double & entryY, 
                                   const double & entryHrad){
  #ifndef NDEBUG
  // check if heading is a valid value
  if ( !validAngle(entryHrad) ){
    char msg[150]; 
    sprintf(msg, "Error, Target::set_entryState: entryHrad invalid: must be in radians on the interval [0,2*pi], but got : %i", entryHrad);
    throw std::runtime_error(msg);
  }
  #endif
  // determine entry angle and direction corresponding to supplied state
  double entryAngleRad = polarAngleFromState(entryX, entryY);
  // check if heading is tangent and determine orbitDir
  int orbitDir;  
  tangentHeading(entryHrad, entryAngleRad, orbitDir);
  // check if the supplied state lies on the orbit
  pointOnOrbit(entryX, entryY);
  // set private variables
  entryX_ = entryX;
  entryY_ = entryY;
  entryHrad_ = entryHrad;
  entryAngleRad_ = entryAngleRad;
  orbitDir_ = orbitDir;
  // determine the exit state
  exitAngleRad_ = MathTools::mod( entryAngleRad_ 
                                + (double)orbitDir_*orbitAngleRad_, 2.0*M_PI);
  exitX_ = centerX_ + orbitRadius_*cos( exitAngleRad_ );
  exitY_ = centerY_ + orbitRadius_*sin( exitAngleRad_ );
  exitHrad_ = MathTools::mod( exitAngleRad_ 
                            + (double)orbitDir*M_PI/2.0, 2.0*M_PI);  
  targetDefined_ = true;
  return true;
}

bool ODTSP::Target::set_exitState(const double & exitX, 
                                  const double & exitY, 
                                  const double & exitHrad){
  #ifndef NDEBUG
  // check if heading is a valid value
  if ( !validAngle(exitHrad) ){
    char msg[150]; 
    sprintf(msg, "Error, Target::set_exitState: exitHrad invalid: must be in radians on the interval [0,2*pi], but got : %i", exitHrad);
    throw std::runtime_error(msg);
  }
  #endif
  // determine exit angle and direction corresponding to supplied state
  double exitAngleRad = polarAngleFromState(exitX, exitY);
  // check if heading is tangent and determine orbitDir
  int orbitDir;  
  tangentHeading(exitHrad, exitAngleRad, orbitDir);
  // check if the supplied state lies on the orbit
  pointOnOrbit(exitX, exitY);
  // set private variables
  exitX_ = exitX;
  exitY_ = exitY;
  exitHrad_ = exitHrad;
  exitAngleRad_ = exitAngleRad;
  orbitDir_ = orbitDir;
  // determine the entry state
  entryAngleRad_ = MathTools::mod( exitAngleRad_  
                                - (double)orbitDir_*orbitAngleRad_, 2.0*M_PI);
  entryX_ = centerX_ + orbitRadius_*cos( entryAngleRad_ );
  entryY_ = centerY_ + orbitRadius_*sin( entryAngleRad_ );
  entryHrad_ = MathTools::mod( entryAngleRad_ 
                            + (double)orbitDir*M_PI/2.0, 2.0*M_PI);  
  targetDefined_ = true;
  return true;
}

bool ODTSP::Target::pointOnOrbit(const double & x,const double & y){
  // distance of point from center squared
  double dsqr = (x-centerX_)*(x-centerX_) + (y-centerY_)*(y-centerY_);  
  // compare to orbit radius squared
  if ( std::abs( dsqr - orbitRadius_*orbitRadius_ ) > tolLength_ ){
    #ifndef NDEBUG
    char msg[150]; 
    sprintf(msg, "Error, Target::pointOnOrbit: entryX and entryY invalid: does not match radius orbitRadius specified. point on orbit (x,y) = (%3.3f, %3.3f), orbit center (x,y) = (%3.3f, %3.3f)", 
            x, y, centerX_, centerY_);
    throw std::runtime_error(msg);
    #endif
    return false;
  }
  else {
    return true;
  }
}

bool ODTSP::Target::tangentHeading(const double & headingRad, 
                                   const double & polarAngleRad, 
                                   int & orbitDir){
  // consider two possible directions of travel
  // the sign of the entry angle determines if the heading corresponds to a 
  // clockwise (CW) or counter-clockwise (CCW) orbit
  double headingRad_CCW = polarAngleRad + M_PI/2.0;
  double headingRad_CW = polarAngleRad - M_PI/2.0;
  // check if this is consistent with the headingRad supplied
  if ( MathTools::polarDistance(headingRad_CCW, headingRad) < tolAngleRad_ ){
    orbitDir = 1; 
    return true;
  }  
  else if ( MathTools::polarDistance(headingRad_CW, headingRad) 
            < tolAngleRad_ ){
    orbitDir = -1;
    return true;
  }  
  else { 
    #ifndef NDEBUG
    char msg[150]; 
    sprintf(msg, "Error, Target::tangentHeading: entryHrad invalid angle is not tangent to the orbit, headingRad: %3.3f, polarAngleRad: %3.3f , headingRad_CCW: %3.3f, headingRad_CW: %3.3f", 
            headingRad, polarAngleRad, headingRad_CCW, headingRad_CW );
    throw std::runtime_error(msg); 
    #endif
    return false;
  }
}

bool ODTSP::Target::validAngle(const double & h){
  if ( h < 0 || h> 2.0*M_PI ){    
    return false;
  }
  else {
    return true;
  }
}

double ODTSP::Target::polarAngleFromState(const double & x, const double & y){
  return MathTools::mod( atan2( y - centerY_ , x - centerX_ ), 2.0*M_PI);  
}

void ODTSP::Target::print(){
  if ( elems_.size() > 0 ){
    if ( type_ == ODTSP::TargetType::Single ){
      printf("Single Element No.: %d", elems_.at(0) );
    }
    else {
      printf("Cluster Elements No.:");
      for (int i = 0; i < elems_.size() ; i++){
        printf("%i ",elems_.at(i));
      }
    }
    printf("\n");
    printf("Center (xc,yc) : (%3.3f, %3.3f)\n",centerX_, centerY_);
    printf("Radius : %3.3f \n", orbitRadius_);
    printf("Orbit Angle (deg) : %3.3f \n", orbitAngleRad_*180/M_PI);
    if ( targetDefined_ ){
      printf("Entry State (x,y,hdeg) = (%3.3f, %3.3f, %3.3f)\n",
                                                               entryX_, entryY_, 
                                                            entryHrad_*180/M_PI); 
      printf("Exit State (x,y,hdeg) = (%3.3f, %3.3f, %3.3f)\n", exitX_, exitY_, 
                                                             exitHrad_*180/M_PI);
      std::string orbitDirString = "CCW";
      if ( orbitDir_ == -1){
        orbitDirString = "CW";
      }
      printf("(entryAngleDeg, exitAngleDeg, orbitDir) = (%3.3f, %3.3f, %i (%s))\n", 
                                entryAngleRad_*180/M_PI ,exitAngleRad_*180/M_PI, 
                                orbitDir_, orbitDirString.c_str() ); 
    }
  }
  else {
    printf("Warning: Target has no data!\n");
  }
}

// get functions: the target must be defined for these to work 
double ODTSP::Target::get_entryX(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_entryX invalid because target is not fully defined. "); 
  }
  #endif
  return entryX_;
};

double ODTSP::Target::get_entryY(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_entryY invalid because target is not fully defined. "); 
  }
  #endif
  return entryY_;
};

double ODTSP::Target::get_entryHrad(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_entryHrad invalid because target is not fully defined. "); 
  }
  #endif
  return entryHrad_;
}; 

double ODTSP::Target::get_entryAngleRad(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_entryAngleRad invalid because target is not fully defined. "); 
  }
  #endif
  return entryAngleRad_;
}; 

double ODTSP::Target::get_exitX(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_exitX invalid because target is not fully defined. "); 
  }
  #endif
  return exitX_;
}; 

double ODTSP::Target::get_exitY(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_exitY invalid because target is not fully defined. "); 
  }
  #endif
  return exitY_;
}; 

double ODTSP::Target::get_exitHrad(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_exitHrad invalid because target is not fully defined. "); 
  }
  #endif
  return exitHrad_;
}; 

double ODTSP::Target::get_exitAngleRad(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_exitAngleRad invalid because target is not fully defined. "); 
  }
  #endif
  return exitAngleRad_;
}; 

double ODTSP::Target::get_centerX(){
  return centerX_;
}; 

double ODTSP::Target::get_centerY(){
  return centerY_;
}; 

double ODTSP::Target::get_orbitRadius(){
  return orbitRadius_;
}; 

double ODTSP::Target::get_orbitAngleRad(){
  return orbitAngleRad_;
}; 

bool ODTSP::Target::contains_elem(int elemInd){
  // check if elemInd is in elems_
  return std::find(elems_.begin(), elems_.end(), elemInd) != elems_.end();
}

bool ODTSP::Target::contains_elem(int elemInd, int & entryInd){
  // check if elemInd is in elems_
  auto pos = std::find(elems_.begin(), elems_.end(), elemInd);
  entryInd = (int)std::distance( elems_.begin(), pos);
  return  pos != elems_.end();
}

int ODTSP::Target::get_orbitDir(){
  #ifndef NDEBUG
  if (!targetDefined_){
    print();
    throw std::runtime_error("Error, Target::get_orbitDir invalid because target is not fully defined. "); 
  }
  #endif
  return orbitDir_;
}; 


vi ODTSP::Target::get_elems(){
  return elems_;
}; 

int ODTSP::Target::get_numElems(){
  return (int)elems_.size();
}

ODTSP::TargetType ODTSP::Target::get_type(){
  return type_;
}; 
