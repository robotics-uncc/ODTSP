#include "ConvectedDubins_Solver.h"
#include "ConvectedDubins_Utils.h"
#include "MathTools.h"
#include<cmath>
#include<algorithm>
#include<limits>

void ConvectedDubins::Solver::solveBSB( ConvectedDubins::pathType pt ){
  double d1, d2;
  if ( pt == ConvectedDubins::pathType::LSL ){
    d1 = -1; // negative delta is a left turn
    d2 = -1;
  }
  else if ( pt == ConvectedDubins::pathType::RSR ){
    d1 = 1;
    d2 = 1;
  }
  else if ( pt == ConvectedDubins::pathType::LSR ){
    d1 = -1; 
    d2 = 1;
  }
  else if ( pt == ConvectedDubins::pathType::RSL ){
    d1 = 1;
    d2 = -1;
  }
  else {
    throw std::runtime_error("ConvectedDubins::Solver::solveBSB: inprob_->v()lid path type.");
  };
  bool solnFound = false;

  // p. 1739, above eq. 20
  double phit1 = fmod( prob_->hInitialRad() - prob_->hwRad(), 2.0*M_PI);
  double phit2 = fmod( prob_->hFinalRad() - prob_->hwRad() - d2*prob_->omega()*t2pi_, 2.0*M_PI);
  // eqs. 20 -23
  double xt10 = x0_ - prob_->v()/(d1*prob_->omega())*sin(phit1);
  double yt10 = y0_ + prob_->v()/(d1*prob_->omega())*cos(phit1);
  double xt20 = x1_ - prob_->v()/(d2*prob_->omega())*sin(d2*prob_->omega()*t2pi_ + phit2)
                                - prob_->w()*t2pi_;
  double yt20 = y1_ + prob_->v()/(d2*prob_->omega())*cos(d2*prob_->omega()*t2pi_ + phit2);
  // analytical solution 
  if ( pt == LSL || pt == RSR ){ 
    for (int kint = -3; kint < 3; kint++){
      double k = (double)kint;
      double num = yt20 - yt10;
      double denom = xt20 - xt10 + prob_->w()*( fmod(phit1-phit2, 2.0*M_PI) - 2.0*k*M_PI )/
                                              ( d1*prob_->omega() );
      double alpha = atan2(num, denom);  
      double tA = t2pi_/(d1*2.0*M_PI)*( asin(prob_->w()/prob_->v()*sin(alpha)) 
                                        + alpha - phit1 ); // ignore m? eq. 36
      if ( tA > t2pi_ || tA < 0){
        tA = tA - t2pi_*floor(tA/t2pi_);
      }
      // eq. 34
      double tB = tA + ( fmod(phit1 - phit2, 2.0*M_PI) - 2.0*k*M_PI ) 
                         / (d1*2.0*M_PI) * t2pi_; // A.W. 2-Jan-2019, t2pi term not in eq. 34?
      if ( checkConditionsBSB( d1, d2, tA, tB, xt10, yt10, phit1, xt20, yt20, phit2, pt ) ){
       // printf("conditions met.\n");
        solnFound = true;
      }
    }
  }
  // numerical solution 
  else { 
    for (int kint = -2; kint < 3; kint++){
      double k = (double)kint;
      // test different initial conditions
      std::vector<double> rootsComputed(numTestPts_);
      for ( int l = 0; l < numTestPts_; l++){
        bool rootAlreadyfound = false;
        // initial guess
        double p0 = (double)(l+1) * t2pi_ / (double)numTestPts_;
        // root solving 
        double tA = fixedPointBSB( p0, d1, d2, k,
                                   xt10, yt10, phit1, 
                                   xt20, yt20, phit2 );
        // check bounds
        if ( tA > t2pi_ || tA < 0){
          tA = tA - t2pi_*floor(tA/t2pi_);
        }
        // check if the root has already been computed 
        if ( l > 0 ){
          for (int i = 0; i < l; i++){
            if ( fabs( rootsComputed[i] - tA ) <= sameRootEpsilon_ ){
              rootAlreadyfound = true;
            } 
          }
        }
        // store the computed root
        rootsComputed[l] = tA;
        // if the root is unique
        if ( !rootAlreadyfound ){
          double tB = d1/d2*tA + (phit1 - phit2 + 2.0*k*M_PI)/(d2*prob_->omega());
          if ( checkConditionsBSB( d1, d2, tA, tB, xt10, yt10, phit1, xt20, yt20, phit2, pt ) ){
            //printf("conditions met.\n");
            solnFound = true;
          }
        }
      }
    }
  }
  if ( !solnFound ){
    switch ( pt ){
      case pathType::LSL: LSL_.set_pathStatusInfeasible(); break;
      case pathType::RSR: RSR_.set_pathStatusInfeasible(); break;
      case pathType::LSR: LSR_.set_pathStatusInfeasible(); break;
      case pathType::RSL: RSL_.set_pathStatusInfeasible(); break;
    }
  }
}

void ConvectedDubins::Solver::solveBBB( ConvectedDubins::pathType pt ){
  double d1, d2;
  if ( pt == ConvectedDubins::pathType::LRL ){
    d1 = -1; // negative delta is a left turn
  }
  else if ( pt == ConvectedDubins::pathType::RLR ){
    d1 = 1;
  }
  else {
    throw std::runtime_error("ConvectedDubins::Solver::solveBBB: invalid path type.");
  };
  d2 = -d1;

  bool solnFound = false;

  // p. 1739, above eq. 20
  double phit1 = fmod( prob_->hInitialRad() , 2.0*M_PI);
  double xt10 = x0_ - prob_->v()/(d1*prob_->omega())*sin(phit1);
  double yt10 = y0_ + prob_->v()/(d1*prob_->omega())*cos(phit1);

  // test a grid of initial conditions, grid resolution
  std::vector< std::vector<double> > rootsComputed(numTestPts_*numTestPts_);
  for (int i = 0; i < rootsComputed.size(); i++){
    rootsComputed[i].resize(3);
  }
  for (int kint = -2; kint < 3; kint++){
    double k = (double)kint;
    for ( int l = 0; l < numTestPts_; l++){
      for ( int n = 0; n < numTestPts_; n++){
        bool rootAlreadyfound = false;
        // initial guess
        std::vector<double> pvecInit(2);
        std::vector<double> pvecOut(2);
        pvecInit[0] = (double)(l+1)*t2pi_/numTestPts_;
        pvecInit[1] = (double)(n+1)*3*t2pi_/numTestPts_;
        fixedPointBBB( pvecInit, d1, k, pvecOut, xt10, yt10, phit1);
        double tA = pvecOut[0];
        double T = pvecOut[1];
        double tB = tA + T/2 + (prob_->hFinalRad() - prob_->hInitialRad() + 2.0*k*M_PI)
                              /( 2*d2*prob_->omega() );
        int curRow = (l)*numTestPts_+(n);
        //printf("testing...root eps_ilon \n");
        // check if the root has already been computed 
        if ( l > 0 || n > 0 ){
          for (int i = 0; i < curRow; i++){
            if (    fabs( rootsComputed[i][0] - tA ) <= sameRootEpsilon_ 
                 && fabs( rootsComputed[i][1] - tB ) <= sameRootEpsilon_ 
                 && fabs( rootsComputed[i][2] - T ) <= sameRootEpsilon_ ){
              rootAlreadyfound = true;
            } 
          }
        }
        // store the computed rootz
        rootsComputed[curRow][0] = tA;
        rootsComputed[curRow][1] = tB;
        rootsComputed[curRow][2] = T;

        // if the root is unique
        if ( !rootAlreadyfound ){
          if ( checkConditionsBBB( d1, tA, tB, T, xt10, yt10, phit1, pt ) ){
            //printf("conditions met.\n");
            solnFound = true;
          }
        }
      }
    }
  }
 
  if ( !solnFound ){
    switch ( pt ){
      case pathType::LRL: LRL_.set_pathStatusInfeasible(); break;
      case pathType::RLR: RLR_.set_pathStatusInfeasible(); break;
    }
  }
  
}

double ConvectedDubins::Solver::fixedPointBSB( double p0, double d1, double d2, double k, 
                                               double xt10, double yt10, double phit1,
                                               double xt20, double yt20, double phit2 ){
  std::vector<double> pvec(N_);
  pvec[0] = p0;
  int i = 1;
  while (i <= N_){
    double p = pvec[i-1];
    double E = (d1-d2) / ( d2*d1*prob_->omega() ) - (yt20-yt10)/prob_->w();
    double F = (xt20-xt10)/prob_->w() + (d1/d2-1)*p + (phit1-phit2 + 2.0*k*M_PI) 
                                                          /(d2*prob_->omega());
    double G = 1/prob_->v()*(yt20 - yt10) + prob_->v()/prob_->w()*(d2-d1)/(d1*d2*prob_->omega());
    double f = E*cos(d1*prob_->omega()*p + phit1) + F*sin(d1*prob_->omega()*p + phit1) - G;
    double fBar = - E*sin(d1*prob_->omega()*p + phit1)*d1*prob_->omega() 
                  + F*cos(d1*prob_->omega()*p + phit1)*d1*prob_->omega()
                  + sin(d1*prob_->omega()*p + phit1)*(d1/d2-1);
    pvec[i] = p - f/fBar;
    if ( fabs( pvec[i] - pvec[i-1]) < tol_ ){
      return pvec[i];
    }
    i++;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

bool ConvectedDubins::Solver::fixedPointBBB( std::vector<double> pvecInit, double d1, double k, 
                                             std::vector<double> & pvecOut,
                                             double xt10, double yt10, double phit1){
  pvecOut.resize(2);
  double d2 = -d1;
  double d3 = d1;
  std::vector<double> pvec1(N_);
  std::vector<double> pvec2(N_);
  pvec1[0] = pvecInit[0];
  pvec2[0] = pvecInit[1];
  int i = 1;
  while (i <= N_){
    double tA = pvec1[i-1];
    double T = pvec2[i-1];
    // vector F = [f1; f2];
    double f1 = 2*prob_->v()/(d1*prob_->omega())*sin(d1*prob_->omega()*tA+phit1) 
                + prob_->w()*T +xt10 - prob_->xFinal() + prob_->v()/(d3*prob_->omega())
                *sin(prob_->hFinalRad()) + 2*prob_->v()/(d2*prob_->omega())
                *sin( d2*prob_->omega()*T/2 + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                      + d1*prob_->omega()*tA );
    double f2 = -2*prob_->v()/(d1*prob_->omega())*cos(d1*prob_->omega()*tA+phit1) 
                  + yt10 - prob_->yFinal() - prob_->v()/(d3*prob_->omega())*cos(prob_->hFinalRad()) 
                  - 2*prob_->v()/(d2*prob_->omega())
                  *cos( d2*prob_->omega()*T/2 + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                        + d1*prob_->omega()*tA );
    // matrix FBar = [ a b ; c d ];
    double a = 2*prob_->v()*( cos(d1*prob_->omega()*tA + phit1) 
                             - cos( d2*prob_->omega()*T/2 + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                                    + d1*prob_->omega()*tA ) );
    double b = prob_->w()+prob_->v()*cos( d2*prob_->omega()*T/2 
                                          + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                                          + d1*prob_->omega()*tA );
    double c = 2*prob_->v()*( sin(d1*prob_->omega()*tA + phit1) 
                             - sin( d2*prob_->omega()*T/2 + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                                    + d1*prob_->omega()*tA ) );
    double d = prob_->v()*sin( d2*prob_->omega()*T/2 + (prob_->hFinalRad()+phit1 + k*2*M_PI)/2 
                               + d1*prob_->omega()*tA );
    double det = a*d - b*c;
    // matrix FBarInv = [d/det -b/det; -c/det a/det];
    // pvec = x - FBar*FBarInv;
    pvec1[i] = pvec1[i-1] - (d/det*f1 - b/det*f2);
    pvec2[i] = pvec2[i-1] - (-c/det*f1 + a/det*f2);
    // convergence criteria
    if ( (pvec1[i]-pvec1[i-1])*(pvec1[i]-pvec1[i-1]) 
         + (pvec2[i]-pvec2[i-1])*(pvec2[i]-pvec2[i-1]) < tol_ ){
      pvecOut[0] = pvec1[i];
      pvecOut[1] = pvec2[i];
      return true;
    }
  i++;
  }
  pvecOut[0] = std::numeric_limits<double>::quiet_NaN();
  pvecOut[1] = std::numeric_limits<double>::quiet_NaN();
  return false;
}


bool ConvectedDubins::Solver::straightLineSolnExists(){

  // determine if final heading is equal to initial heading within tolerance 
  double headingDiffRad = MathTools::polarDistance( prob_->hInitialRad() , prob_->hFinalRad() );
  if ( headingDiffRad <= straightLineCourseTolRad_ ){ 
    // compute current course angle 
    double initialCourseRad = courseAngleRad( prob_->v() , prob_->w() , prob_->hInitialRad() );
    // compute angle to desired point 
    double delx = prob_->xFinal() - prob_->xInitial();
    double dely = prob_->yFinal() - prob_->yInitial();
    double courseToEndpointRad = atan2( dely, delx );

    // determine if point is along velcotiy vector within tolerances    
    double courseDiffRad = MathTools::polarDistance( initialCourseRad , courseToEndpointRad );
    if ( courseDiffRad <= straightLineCourseTolRad_ ){ 
      // make LSL the solution 
      LSL_.set_prob( prob_ );
      LSL_.set_pathType( LSL );
      ConvectedDubins::Path * path; 
      path = &LSL_;
      // speed over ground, first decompose 
      double vNorth = prob_->v() * cos( initialCourseRad ) + prob_->wx();
      double vEast = prob_->v() * sin( initialCourseRad ) + prob_->wy();
      double vg = sqrt(vNorth*vNorth + vEast*vEast);
      double pathLength = sqrt(delx*delx + dely*dely);
      double T  = pathLength / vg; // time to traverse 
      path->set_pathParamsBSB(1, 1, 0.0, T, T); 
      return true;
    }
  }
  return false;
}

bool ConvectedDubins::Solver::checkConditionsBSB( double d1, double d2, double tA, double tB, 
                                                  double xt10, double yt10, double phit1,
                                                  double xt20, double yt20, double phit2,
                                                  ConvectedDubins::pathType pt ){
  // check that tA, tB, and T are valid numbers 
  if ( std::isnan(tA) || std::isnan(tB) ){
    return false;
  }

  // obtained by differentiating eq. 18,19
  double xtAdot = prob_->v()*cos(d1*prob_->omega()*tA + phit1) + prob_->w();
  double ytAdot = prob_->v()*sin(d1*prob_->omega()*tA + phit1);
  double xtBdot = prob_->v()*cos(d2*prob_->omega()*tB + phit2) + prob_->w();
  double ytBdot = prob_->v()*sin(d2*prob_->omega()*tB + phit2);
  // eq. 16
  double xtA = prob_->v()/(d1*prob_->omega())*sin(d1*prob_->omega()*tA + phit1)
                            + prob_->w()*tA + xt10;
  // eq. 17
  double ytA = -prob_->v()/(d1*prob_->omega())*cos(d1*prob_->omega()*tA + phit1) + yt10;
  // eq. 18
  double xtB = prob_->v()/(d2*prob_->omega()) *sin(d2*prob_->omega()*tB + phit2)
                            + prob_->w()*tB + xt20;
  // eq. 19
  double ytB = -prob_->v()/(d2*prob_->omega())*cos(d2*prob_->omega()*tB + phit2) + yt20;
  // check that range of tB is valid
  if ( tB < -t2pi_ || tB > t2pi_ ){
    return false;
  } 
  // check that range of tA is valid 
  if ( tA < 0 || tA > 2*t2pi_ ){
    //printf(" tA < 0 || tA > 2*t2pi_ \n");
    return false;
  }
  // check that not negative 
  if ( (xtB-xtA)*xtAdot < -eps_ ){
    //printf(" (xtB-xtA)*xtAdot < -eps_ \n");
    return false;
  }  
  if ( (ytB-ytA)*ytAdot < -eps_ ){
    //printf(" (ytB-ytA)*ytAdot < -eps_ \n"); 
    return false;
  } 
  if ( (xtB-xtA)*xtBdot < -eps_ ){
    //printf(" (xtB-xtA)*xtBdot < -eps_ \n"); 
    return false;
  }  
  if ( (ytB-ytA)*ytBdot < -eps_ ){
    //printf(" (ytB-ytA)*ytBdot < -eps_ \n"); 
    return false;
  } 
  // check that equal 
  if ( fabs( (xtB - xtA)/(ytB-ytA) - xtAdot/ytAdot ) >= eps_ ){
    //printf(" fabs( (xtB - xtA)/(ytB-ytA) - xtAdot/ytAdot ) >= eps_ \n"); 
    return false;
  }  
  if ( fabs( (xtB - xtA)/(ytB-ytA) - xtBdot/ytBdot ) >= eps_ ){
    //printf("fabs( (xtB - xtA)/(ytB-ytA) - xtBdot/ytBdot ) >= eps_ \n"); 
    return false;
  }  
  // check endpoint
  // from prop. 1
  double tBeta = tA + sqrt( (xtB-xtA)*(xtB-xtA) + 
                            (ytB-ytA)*(ytB-ytA) )
                     /sqrt( (xtBdot*xtBdot + ytBdot*ytBdot) ); 
  double T = tBeta + (t2pi_ - tB);  

  // check if the candidate satisfies the desired final (x,y) position 
  double xFinal, yFinal, hFinalRad;
  computeEndpoint( prob_, pt, d1, d2, tA, tBeta, T, xFinal, yFinal, hFinalRad);
  if (    fabs(xFinal - prob_->xFinal()) > eps_ 
       && fabs(yFinal - prob_->yFinal()) > eps_ ){
    //printf("endpoint not equal. \n"); 
    return false;
  }
  
  ConvectedDubins::Path * path; 
  switch ( pt ){
    case pathType::LSL: path = &LSL_; break;
    case pathType::RSR: path = &RSR_; break;
    case pathType::LSR: path = &LSR_; break;
    case pathType::RSL: path = &RSL_; break;
  }
  // assign the params if the path is not yet defined, or if it is defined and 
  // the old T is larger than the new T
  if ( !path->isDefined() || ( path->isDefined() && path->get_T() > T) ){
    path->set_pathParamsBSB(d1, d2, tA, tBeta, T); 
  } 
  return true;
}

bool ConvectedDubins::Solver::checkConditionsBBB( double d1, double tA, double tB, double T,
                                                  double xt10, double yt10, double phit1,
                                                  ConvectedDubins::pathType pt ){

  // check that tA, tB, and T are valid numbers 
  if ( std::isnan(tA) || std::isnan(tB) || std::isnan(T) ){
    return false;
  }

  double d2 = -d1;
  double d3 = d1;

  double phit3 = fmod ( prob_->hFinalRad() - d3*prob_->omega()*T , 2.0*M_PI);
  double xt30 = prob_->xFinal() - prob_->v()/(d3*prob_->omega())*sin( prob_->hFinalRad() ) - prob_->w()*T;
  double yt30 = prob_->yFinal() + prob_->v()/(d3*prob_->omega())*cos( prob_->hFinalRad() );

  double phit2 = 2*d1*prob_->omega()*tA + phit1;
  double xt20 = xt30 - 2*prob_->v()/(d2*prob_->omega())*sin(d2*prob_->omega()*tB + phit2);
  double yt20 = yt30 + 2*prob_->v()/(d2*prob_->omega())*cos(d2*prob_->omega()*tB + phit2);


  // eq. 16-19
  double xt1tA = prob_->v()/(d1*prob_->omega())*sin(d1*prob_->omega()*tA + phit1)
                            + prob_->w()*tA + xt10;
  double yt1tA = -prob_->v()/(d1*prob_->omega())*cos(d1*prob_->omega()*tA + phit1) + yt10;

  double xt2tA = prob_->v()/(d2*prob_->omega())*sin(d2*prob_->omega()*tA + phit2)
                            + prob_->w()*tA + xt20;
  double yt2tA = -prob_->v()/(d2*prob_->omega())*cos(d2*prob_->omega()*tA + phit2) + yt20;

  //
  double xt2tB = prob_->v()/(d2*prob_->omega())*sin(d2*prob_->omega()*tB + phit2)
                            + prob_->w()*tB + xt20;
  double yt2tB = -prob_->v()/(d2*prob_->omega())*cos(d2*prob_->omega()*tB + phit2) + yt20;
  double xt3tB = prob_->v()/(d3*prob_->omega())*sin(d3*prob_->omega()*tB + phit3)
                            + prob_->w()*tB + xt30;
  double yt3tB = -prob_->v()/(d3*prob_->omega())*cos(d3*prob_->omega()*tB + phit3) + yt30;  

  // derived by hand by differentiating 18,19
  double xt1dottA = prob_->v()*cos(d1*prob_->omega()*tA + phit1) + prob_->w();
  double yt1dottA = prob_->v()*sin(d1*prob_->omega()*tA + phit1);
  double xt2dottA = prob_->v()*cos(d2*prob_->omega()*tA + phit2) + prob_->w();
  double yt2dottA = prob_->v()*sin(d2*prob_->omega()*tA + phit2);

  double xt2dottB= prob_->v()*cos(d2*prob_->omega()*tB + phit2) + prob_->w();
  double yt2dottB = prob_->v()*sin(d2*prob_->omega()*tB + phit2);
  double xt3dottB = prob_->v()*cos(d3*prob_->omega()*tB + phit3) + prob_->w();
  double yt3dottB = prob_->v()*sin(d3*prob_->omega()*tB + phit3);

  // check that: 0 <= tA <= tB <= T <= 4.0*t2pi   
  if ( tA < 0 || tB < tA || T < tB || T > 4.0*t2pi_ ){
    //printf(" tA < 0 || tB < tA || T < tB || T > 4.0*t2pi_ \n");
    return false;
  }
  if ( fabs(xt1tA - xt2tA) >= eps_ ){
    //printf(" fabs(xt1tA - xt2tA) >= eps_ \n");
    return false; 
  }
  if ( fabs(yt1tA - yt2tA) >= eps_ ){
    //printf("  fabs(yt1tA - yt2tA) >= eps_ \n");
    return false; 
  }
  if ( fabs(xt2tB - xt3tB) >= eps_ ){
    //printf(" fabs(xt2tB - xt3tB) >= eps_ \n");
    return false; 
  }
  if ( fabs(yt2tB - yt3tB) >= eps_ ){
    //printf(" fabs(xt2tB - xt3tB) >= eps_ \n");
    return false; 
  }
  if ( fabs(xt1dottA - xt2dottA) >= eps_ ){
    //printf(" fabs(xt1dottA - xt2dottA) >= eps_ \n");
    return false; 
  }
  if ( fabs(yt1dottA - yt2dottA) >= eps_ ){
    //printf(" fabs(yt1dottA - yt2dottA) >= eps_ \n");
    return false; 
  }
  if ( fabs(xt2dottB - xt3dottB) >= eps_ ){
    //printf(" fabs(xt2dottB - xt3dottB) >= eps_ \n");
    return false; 
  }
  if ( fabs(yt2dottB - yt3dottB) >= eps_ ){
    //printf(" fabs(yt2dottB - yt3dottB) >= eps_ \n");
    return false; 
  }
  // check if the candidate satisfies the desired final (x,y) position 
  double xFinal, yFinal, hFinalRad;
  computeEndpoint( prob_, pt, d1, d3, tA, tB, T, xFinal, yFinal, hFinalRad);
  if (    fabs(xFinal - prob_->xFinal()) > eps_ 
       && fabs(yFinal - prob_->yFinal()) > eps_ ){
    return false;
  }
  ConvectedDubins::Path * path; 
  switch ( pt ){
    case pathType::LRL: path = &LRL_; break;
    case pathType::RLR: path = &RLR_; break;
  }
  // assign the params if the path is not yet defined, or if it is defined and the old T is larger
  // than the new T
  if ( !path->isDefined() || ( path->isDefined() && path->get_T() > T) ){
    path->set_pathParamsBBB(d1, tA, tB, T);
  } 
  return true;
}

ConvectedDubins::Solver::Solver( ConvectedDubins::Problem * prob ){
  #ifndef NDEBUG
  if ( !prob->isDefined() ){
    throw std::runtime_error("ConvectedDubins::Solver::Solver: Problem not defined.");
  }
  #endif 
  prob_ = prob;
  // Lemma 2, p. 1738
  t2pi_ = 2.0*M_PI/prob_->omega();
  ConvectedDubins::convertInertialtoTrochoid( prob_->xInitial(), 
                                              prob_->yInitial(),
                                              prob_->hwRad(),
                                              x0_, y0_ );
  ConvectedDubins::convertInertialtoTrochoid( prob_->xFinal(), prob_->yFinal(),
                                              prob_->hwRad(),
                                              x1_, y1_ );
  // check speciale case
  if ( !straightLineSolnExists() ){
    LSL_.set_prob( prob_ );
    LSL_.set_pathType( LSL );
    
    RSR_.set_prob( prob_ );
    RSR_.set_pathType( RSR );
    
    RSL_.set_prob( prob_ );
    RSL_.set_pathType( RSL );
    
    LSR_.set_prob( prob_ );
    LSR_.set_pathType( LSR );

    RLR_.set_prob( prob_ );
    RLR_.set_pathType( RLR );

    LRL_.set_prob( prob_ ); 
    LRL_.set_pathType( LRL );

    solveBSB( LSL );
    solveBSB( RSR );

    solveBSB( RSL );
    solveBSB( LSR );
   
    solveBBB( RLR );
    solveBBB( LRL );    
  }
  setOptimal();
}

void ConvectedDubins::Solver::setOptimal(){
  // Tvec is a vector of path length (time) for each candidate, intiialized to double max
  std::vector<double> Tvec(6);
  std::fill ( Tvec.begin(), Tvec.end(), std::numeric_limits<double>::max() );
  // if the path is defined, redefine entry from double max to the computed T value
  if ( LSL_.isDefined() ){
    Tvec[0] = LSL_.get_T();
  }
  if ( RSR_.isDefined() ){
    Tvec[1] = RSR_.get_T();
  }
  if ( LSR_.isDefined() ){
    Tvec[2] = LSR_.get_T();
  }
  if ( RSL_.isDefined() ){
    Tvec[3] = RSL_.get_T();
  }
  if ( LRL_.isDefined() ){
    Tvec[4] = LRL_.get_T();
  }
  if ( RLR_.isDefined() ){
    Tvec[5] = RLR_.get_T();
  }
  // determine which T value is the smallest 
  auto minIndexIter = std::min_element( std::begin(Tvec) , std::end(Tvec) );
  int minIndex = std::distance( Tvec.begin() , minIndexIter );
  // note the pathType enumeration is defined such that LSL = 0, RSR = 1, etc. corresponding to 
  // the sequence in which Tvec is defined so the following switch case can be used
  // here we point optPath_ to the lowest cost path.
  switch ( minIndex ){
    case pathType::LSL: optPath_ = &LSL_; break;
    case pathType::RSR: optPath_ = &RSR_; break;
    case pathType::LSR: optPath_ = &LSR_; break;
    case pathType::RSL: optPath_ = &RSL_; break;
    case pathType::LRL: optPath_ = &LRL_; break;
    case pathType::RLR: optPath_ = &RLR_; break;
  }
  optPath_->set_pathStatusOptimal();
}




ConvectedDubins::Path ConvectedDubins::Solver::get_path( ConvectedDubins::pathType pt ){
  switch ( pt ){
    case pathType::LSL: return LSL_;
    case pathType::RSR: return RSR_;
    case pathType::LSR: return LSR_;
    case pathType::RSL: return RSL_;
    case pathType::RLR: return RLR_;
    case pathType::LRL: return LRL_;
  }
}

void ConvectedDubins::Solver::computeConnectingStates(){
  LSL_.computeConnectingStates();
  RSR_.computeConnectingStates();
  LSR_.computeConnectingStates();
  RSL_.computeConnectingStates();
  LRL_.computeConnectingStates();
  RLR_.computeConnectingStates();
}

void ConvectedDubins::Solver::plotSolns(std::string fileName){
  LSL_.computePathHistory();
  LSL_.plotXY(fileName,1);
  RSR_.computePathHistory();
  RSR_.plotXY(fileName,1);
  RSL_.computePathHistory();
  RSL_.plotXY(fileName,1);
  LSR_.computePathHistory();
  LSR_.plotXY(fileName,1);
  LRL_.computePathHistory();
  LRL_.plotXY(fileName,1);
  RLR_.computePathHistory();
  RLR_.plotXY(fileName,1);
}

void ConvectedDubins::Solver::printSolns(){
  printf(" =================================================\n");
  printf(" ============ ConvectedDubins::Solver ============\n");
  printf(" =================================================\n");
  LSL_.print();
  RSR_.print();
  RSL_.print();
  LSR_.print();
  RLR_.print();
  LRL_.print();
}

void ConvectedDubins::Solver::get_optimalWaypoints( vd & x, vd & y, vd & hRad ){
  optPath_->get_waypoints(x,y,hRad);
}

void ConvectedDubins::Solver::get_optimalWaypoints_XYEastNorth( vd & x, vd & y, vd & hRad ){
  optPath_->get_waypoints_XYEastNorth(x,y,hRad);
}


