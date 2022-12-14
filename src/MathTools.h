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

// standard headers
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<iomanip>      // std::setprecision
#include<iostream> // std::cout
#include<cmath> // std::sqrt
#include<algorithm> //std::min_element
#include<time.h>
#include<stdlib.h>
#include<random>
// external 
#include<armadillo>

#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H

namespace MathTools {

// PI is given as M_PI
// angles and conversions
const double RAD2DEG = 180.0/M_PI;
const double DEG2RAD = M_PI/180.0;

// return the sign of a number
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// return the modulus of a number, e.g., mod(3pi/2, pi) = pi/2
double        mod(double number, double base); 
int           mod(int number, int base); 
unsigned int  mod(unsigned int number, unsigned int base); 

// determine (smallest) angular angular distance between two angles (rad)
double polarDistance(double a, double b);
double polarDistanceSigned(double a, double b);

double distance(std::vector<double> referenceVector, 
		                       std::vector<double> testVector);

std::vector<int> minIndicesWithTolerance( arma::vec testVector, 
                                          double tolerance);

void append(std::vector<double> &baseVector, std::vector<double> &appendVector);

double sum(std::vector<double> xValues);

// Generates octave commands to plot a curve (with std::vector)
void plot1DCurve(std::string fileName, int figureNumber, 
	               std::vector<double> xData, std::vector<double> yData, 
                 std::string xVarName, std::string yVarName, std::string style);
void plot1DCurve(std::string fileName, int figureNumber, 
                 std::vector<double> xData, std::string xVarNAme, 
                 std::string style);	

// Generates octave commands to plot a curve (with reference to existing vectrs)
void plot1DCurve(std::string fileName, int figureNumber, std::string xVarName, 
	               std::string yVarName, std::string style);
void plot1DCurve(std::string fileName, int figureNumber, std::string xVarName, 
	               std::string style);

void axisProperty(std::string fileName, std::string property);

// Generates octave commands to write vector data to .m file
void writeVectorData(std::string fileName, std::vector<double> data, 
	                   std::string varName);	
void writeVectorData(std::string fileName, std::vector<int> data, 
	                   std::string varName);	
void writeVectorData(std::string fileName, std::vector<double> data, 
	                   std::string varName, int precision);	

// Generate a sequence of evenly spaced points within prescribed bounds
// if desired skip nSkip several points from the begining
std::vector<double> linspace(double xmin, double xmax, double spacing, 
							               int nSkip);
std::vector<double> linspace(double xmin, double xmax, int npts);
std::vector<double> linspace(double xmin, double xmax, double spacing);
std::vector<double> linspace(double xmin, double xmax, int npts, int nSkip);

std::vector<double> polarspace(double angleInitial, double angleFinal, int npts, 
							   int direction);

// Generates (x,y) coordinates defining a circualr arc
	// arcParams is a vector with the following fields:
	// arcParams(0) = xc, x coord of center of arc
	// arcParams(1) = yc, y coord of center of arc
	// arcParams(2) = R, radius of arc
	// arcParams(3) = angleInit, initial angle relative to horizontal
	// arcParams(4) = angleFinal, final angle relative to horizontal
	// numPts = points defining the arc
	// direction = 1, CCW is -1 CW
std::vector< std::vector<double> >  generateCircularArc(
                                            std::vector<double> arcParams, 
                                            int numPts, int direction);

// Generates (x,y) coordinates defining a circle
	// circleParams is a vector with the following fields:
	// circleParams(0) = xc, x coord of center of circle
	// circleParams(1) = yc, y coord of center of circle
	// circleParams(2) = R, radius of circle
std::vector< std::vector<double> > generateCircle(
                                          std::vector<double> circleParams, 
                                          int numPts);

// check if a double value is an integer
bool isInteger(double value);

// Add command
void addCommand(std::string fileName, std::string command);

// Runs a (.m) script in octave
void runOctaveScript(std::string fileName);
void runOctaveScript(std::string fileName, std::string options);

// Returns a random double from a bounded uniform distribution
double randomNumberUniformDistribution(double xmin, double xmax);

// Returns a random vector (std::vector) from a bounded uniform distribution
std::vector<double> randomNumberUniformDistribution(double xmin, double xmax, 
                                                    int npts);

double fmodPos( double val, double mod );


}

#endif
