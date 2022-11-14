#ifndef CONVECTED_DUBINS_PROBLEM_H
#define CONVECTED_DUBINS_PROBLEM_H

#include<vector>

typedef std::vector<double> vd;

namespace ConvectedDubins {

// The ConvectedDubins::Problem object is an input for ConvectedDubins::Solver 
class Problem {

private:

  // Note: all angles are measured in radians, in the NED frame:

  // North-East-Down Reference Frame
  // X : Points North
  // Y : Points East
  // Angles are measured clockwise positive from the X axis 

  // input: vehicle properties
	double R_; // minimum turning radius, default=1
  double v_; // vehicle speed, default=1 
  bool vehPropInputFlag_; // flag if user inputs a non-default radius

  // input: environment properties
  double w_; // wind magnitude, default=0.5 
  double hwRad_; // wind direction, default=PI/2
  double wx_; // wind x-component
  double wy_; // wind y-component
  bool envPropInputFlag_; // flag if user inputs envrionment

	// input: initial condition 
	double xInitial_; // default=0
	double yInitial_; // default=0
	double hInitialRad_; // default=0
  bool startPtInputFlag_; // flag if user inputs a non-default init. cond.

  // input: final condition 
	double xFinal_; // default=0 
	double yFinal_; // default=0 
	double hFinalRad_; // default=0 
  bool endPtInputFlag_; // flag if user inputs endpoint (required)

  // derived quantity
  double omega_; // derived quantity: vehicle maximum turn rate (rad/s)
                 // = v/R

  // called by the constructor 
  void initialize();

  // private set functions
	void set_xInitial( double xInitial );
	void set_yInitial( double yInitial );
	void set_hInitialRad( double hInitialRad );

	void set_xFinal( double xFinal ); 
	void set_yFinal( double yFinal );
	void set_hFinalRad( double hFinalRad );

  // used for displaying startPtInputFlag values, etc. in print()
  void printInputType( bool userInputFlag );

public:

	// constructor and destructor
  Problem( );
	virtual ~Problem(){};

	// set functions
  // note: set these functions before specifying initial/final conditions 
  void set_vehicleProperties( double minTurningRadius, 
                              double nominalSpeed );

  void set_windMagDir( double windMagnitude, double windDirRad);
  void set_windMagXY( double windMagX,  double windMagY);


  // set initial condition 
	void set_stateInitial( double xInitial, 
                         double yInitial, 
                         double hInitialRad );

  // provided as vector of doubles [x y h]
	void set_stateInitial( const vd & stateInitial );

  // set initial conditon in ENU frame ( internally converted to NED frame )
	void set_stateInitial_XYEastNorth( double xInitial, 
                                     double yInitial, 
                                     double hInitialRad );

  // set initial condition in ENU frame with course angle instead of heading
  // This is necessary, for example, if inertial path requires a particular
  // tangent angle. 
	void set_stateInitial_XYEastNorthCourse( double xInitial, 
                                           double yInitial, 
                                           double chiInitialRad );
  
  // see initial condition functions above
	void set_stateFinal( double xFinal, 
                       double yFinal, 
                       double hFinalRad );

  void set_stateFinal( const vd & stateFinal );

	void set_stateFinal_XYEastNorth( double xFinal, 
                                   double yFinal, 
                                   double hFinalRad );

	void set_stateFinal_XYEastNorthCourse( double xFinal, 
                                         double yFinal, 
                                         double chiFinalRad );

	void print(); // prints problem summary to screen 
	
  // there are four flags that must be set for problem to be defined:
  // vehPropInputFlag_, envPropInputFlag_ , startPtInputFlag_, endPtInputFlag_
  bool isDefined();

  // same as above, but omits requirement for endPtInputFlag_
  bool isDefinedNoEndpoint();

  // get functions
	double xInitial(){return xInitial_;}; 
	double yInitial(){return yInitial_;};
	double hInitialRad(){return hInitialRad_;};
	vd stateInitial();

	double xFinal();
	double yFinal();
	double hFinalRad();
	vd stateFinal();

  double R(){return R_;};
  double v(){return v_;};
  double omega(){return omega_;};

  double w(){return w_;};
  double hwRad(){return hwRad_;};
  double wx(){return wx_;};
  double wy(){return wy_;};

}; // class

} // namespace

#endif


