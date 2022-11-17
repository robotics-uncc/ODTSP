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

#ifndef ROBUST_DUBINS_PROBLEM_H
#define ROBUST_DUBINS_PROBLEM_H

#include<iostream> 
#include<vector> 
#include<string> 

#include "RobustDubins_Utils.h"

namespace RobustDubins {

// The RobustDubins::Problem object is mainly used to define a problem as an 
// input for the RobustDubins::Solver 

class Problem {

private:

  // required inputs 
	double m_xFinal; // Final x position
	double m_yFinal; // Final y position
	double m_hFinal; // Final theta (heading) angle [0,2*pi)

	// optional inputs
	double m_xInitial; // default is x0 = 0
	double m_yInitial; // default is y0 = 0
	double m_hInitial; // default is h0 = 0
	double m_minTurningRadius; // default is R=1

  // auto-generated flags
  bool m_minTurningRadiusInputFlag; // flags if user inputs a non-default radius
  bool m_startPtInputFlag; // flags if user inputs a non-default init. cond.
  bool m_endPtDefined; // flag determining if problem is defined 

  // called by the constructor 
  void initialize();

  // we make these set functions private to force user to provide complete
  // initial and final state using set_stateInitial or set_stateFinal
	void set_xInitial( const double & xInitial );
	void set_yInitial( const double & yInitial );
	void set_hInitial( const double & hInitial );

	void set_xFinal( const double & xFinal); 
	void set_yFinal( const double & yFinal);
	void set_hFinal( const double & hFinal);

public:

	// constructor and destructor
  Problem( );
	virtual ~Problem(){};

	// set functions
  // required: (x1,y1,h1) the desired final state defining the problem 
	void set_stateFinal( const double & xFinal, 
                       const double & yFinal, 
                       const double & hFinal );
	void set_stateFinal( const vd & stateFinal );


  // optional: default is (x0,y0,h0) = (0,0,0)
	void set_stateInitial( const double & xInitial, 
                         const double & yInitial, 
                         const double & hInitial );
	void set_stateInitial( const vd & stateInitial );


  // optional: default is R=1
	void set_minTurningRadius( const double & minTurningRadius );


	// main functions
  bool isDefined(); // returns true if the end point has been defined 
	void print(); // prints problem summary to screen 

	// get functions
	double get_xInitial(){return m_xInitial;}; 
	double get_yInitial(){return m_yInitial;};
	double get_hInitial(){return m_hInitial;};
	vd get_stateInitial();

	double get_xFinal();
	double get_yFinal();
	double get_hFinal();
	vd get_stateFinal();

  double get_minTurnRadius(){return m_minTurningRadius;};

  bool get_startPtInputFlag(){return m_startPtInputFlag;};
  bool get_minTurningRadiusInputFlag(){return m_minTurningRadiusInputFlag;};

}; // class

} // namespace

#endif


