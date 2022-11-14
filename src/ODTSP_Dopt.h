#ifndef Dopt_HPP_
#define Dopt_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <memory>
#include <cstring>
#include <stdio.h>
#include <algorithm>
// custom
#include "ODTSP_Problem.h"
#include "ODTSP_Path.h"
#include "ODTSP_Definitions.h"
#include "ODTSP_OrbitList.h"

namespace ODTSP {

enum DoptSolvers{
  ETSP_NO_CLUSTERS,
  GLKH_NO_CLUSTERS,
  GLKH_INTERSECTING_CLUSTERS,
  LKH_NO_CLUSTERS,
  LKH_INTERSECTING_CLUSTERS
};

class Dopt {

  private:
    // vectors indicating tour sequence
    vi etspSoln_; 
    vi atspSoln_; 
    vi gtspSoln_;
    
    // vectors storing states along sthe orbits (positions, headings)
    vd entryAngleRad_;
    vi orbitDir_;

    // specified by user upon construction
    ODTSP::Problem* prob_;  
    ODTSP::OrbitList* list_;

    std::string name_; // name used to initialize

    int ne_; // number of entry states = numEntryAngles*2 
    int nl_; // size of the list
    int nn_; // number of nodes total 

    // applies to intersecting GTSP
    int nt_; // number of single targets
    int nc_; // number of clusters (used in overlapping)
    int nr_; // number of replicas
    vvi Sx_; // the set of clusters which each node belongs to 
    vi targetsInClusters_;

    vi nodeListIndex_;
    vi clusterSize_;
    vi pastReplicas_;

    int precision_;

    // derived from the ODTSP::Problem object (integer variants)
    
    int turnRadius_;
    int orbitRadius_;
    double orbitAngleRad_; // radians

    // intermediate variables
    ODTSP::DoptSolvers sol_;
    int etspCost_; // cost of euclidean tsp
    bool convertedFlag; // flag if converted to ATSP
    int ** C_; // cost matrix
    int M_; // a large number (current max integer / 2)
    int offset_;
    int maxVal_;

    // intermediate functions
    void initialize(); // defines targListIntX/Y, entryAngleRad, orbitDir vars. 

    // use LKH to solve the problem given in _euc.par and saves output to 
    // _euc.tour 
    void solveATSP();

    // populates the C_ with the dubins length cost 
    void buildCostMatrix();

    // modifies the C_ by assigning zero edge costs, or adding a large
    // cost M where appropriate 
    void convertToATSP();

    // write the .gtsp file  
    void writeGTSPfile();

    // write the .tsp and .par file 
    void writeAsymTspParFile();

    // use GLKH to solve the .par file, read in the .tour output and 
    // save it as the the atspSoln_
    void solveGTSP();

    // define gtspSoln_ the same as atspSoln_ but forming a closed loop 
    void extractSoln(); 

    // writes the .sol, _euc.sol, and _euc.cfg files 
    void printSoln();

    // optional outputs

    void resizeSolnVectors( const int & n);

    vd fromX_, fromY_, fromHrad_;
    vi fromCluster_, fromNode_; 

    vd toX_, toY_, toHrad_;
    vi toCluster_, toNode_; 

    int dubinsCost( const int & fromTargInd,
                    const int & fromNodeInd, 
                    const int & toTargInd,
                    const int & toNodeInd );

    int orbitCost(  const int & toTargInd, const int & toNodeInd );

    // check if ATSP solution node i is a replica with node j 
    bool areReplicasSameCluster( const int & iNode, const int & jNode );
    bool isReplica( const int & node );
    int whichCluster( const int & node );
    
    void plotCostMatrixMatlab();

    // ------------------------------------------------------
    void buildCostMatrixIntersecting();
  

    

  public:
    // precision: scale all the length units since we are using integers only
    Dopt(ODTSP::Problem* prob, ODTSP::OrbitList* list, 
         int numPointsPerCluster, int precision, 
         std::string filename);

    // ETSP 
    Dopt(ODTSP::DoptSolvers sol, 
         ODTSP::Problem* prob, ODTSP::OrbitList* list, 
         int precision, std::string filename);
    vi get_etspSoln(){return etspSoln_;};


    ~Dopt();

    // main functions
    void solve(ODTSP::DoptSolvers sol);
    void print();
    void printCostMatrix();
    void printCostMatrixRaw();
    void printGraphProperties();

    // get functions
    
    vi get_gtspSoln(){return gtspSoln_;};
    vi get_atspSoln(){return atspSoln_;};

    vd get_fromX(){return fromX_;};
    vd get_fromY(){return fromY_;};
    vd get_fromHrad(){return fromHrad_;};
    vi get_fromCluster(){return fromCluster_;};

    vd get_toX(){return toX_;};
    vd get_toY(){return toY_;};
    vd get_toHrad(){return toHrad_;};
    vi get_toCluster(){return toCluster_;};

    ODTSP::Path get_path();
    double computeCost();


}; // class


} // namespace





#endif
