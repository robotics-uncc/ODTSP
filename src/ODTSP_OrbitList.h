#ifndef ODTSP_ORBIT_LIST_H
#define ODTSP_ORBIT_LIST_H

#include<vector>
#include "ODTSP_Definitions.h"
#include "ODTSP_Problem.h"
#include "ODTSP_Target.h"


namespace ODTSP {

enum ClusterMode {
  None, // the default target list without clustering 
  Exhaustive, // all clusters, including overlapping subsets are considered
  SimpleCluster, 
  MaximalClusterTuning,
};

class OrbitList {

private:

  ODTSP::Problem* prob_; // pointer to problem definition 
  ODTSP::ClusterMode mode_; // mode used to generate list_

  std::vector<ODTSP::Target> list_; // contains a list of all ODTSP::Targets

  vi seq_; // integers referring to the order of the list_ items in which 
           // the tour proceeds

  // computes an orbit list given a problem and a cluster mode
  void generateTargetList();
  
  // used for Cliquer interface
  // --------------------------------------------------------------------------
  // write a DIMACS Graph File in ASCII format
  // see p. 17-18 of Cliquer documentation 
  //    graphfile : name of the file to be passed to Cliquer
  //    F : a square adjacency matrix with binary entries indicating edges
  void writeDIMACSfile( std::string graphfile, vvi & F, int numVert, int numEdge );

  // 
  void getAdjacencyMatrix(vvi &F, int & numEdges, double eta);

  // request only the maximal cliques for the graph 
  vvi maximalCliques( std::string graphfile );
  
  // request all of the cliques for the graph 
  vvi allCliques( std::string graphfile );

  void assignMaximalCliquesDisjoint(vvi cliques);

  // used by ClusterMode::SimpleCluster
  // --------------------------------------------------------------------------
  // creates a set of bins that stores cluster candidates of different sizes
  void getNeighborSets( vvi & combos, vvi F);
  // delete all repeated entries from 'combos' (according to the equal_seq 
  // criteria)
  void deleteRepeats(vvi & combos);
  // searches through the binnedCombos for valid clusters, rejecting orbit  
  // or truncating and re-testing as appropriate 
  void breakAndAddUniqueClustersToList( vvvi & binnedCombos );
  // returns true of a candidateCluster can be clustered and satisfy the 
  // distance criteria
  bool isValid( vi & candCluster, ODTSP::Target & targ );
  // if a cluster is invalid, this routine removes the most distance point
  // from the cluster center   
  void removeDistantPoint( vi & cluster );

  // used by ClusterMode::Exhaustive
  // --------------------------------------------------------------------------

  void powerSetsOfMaximalCliques(vvi & comboSet);

  // populates the vvi comboSet with possible clusters (and sub-clusters) 
  // based on connectivity
  void getCandidateClusters( vvi & comboSet , vvi F);
  // assigns a center point and radius to each cluster in the orbit list
  void runMiniballAndAddValidToList(const vvi & comboSet);
  // sort the list_ according to the seq given and truncate if necessary
  bool sortAndTruncate(vi seq);

  void binCombosBySize(vvvi &binnedCombos, vvi &comboSet);

  // check if an index is within the limits of the current list 
  bool validIndex(int ind);
  // check if each entry of a sequence is within the limits of the current list 
  bool validSequence(vi seq);

  // when using maximal clusters tuning
  double tuningParam_;

public: 
  // this class is intended to be created as an input into the path planning
  // algorithm with target centers and radii specified when constructed. 
  // the planning algorithm can then determine which clusters or targets to 
  // orbit and then update the properties of each orbit (and rearrange).
  OrbitList( ODTSP::Problem* prob, ODTSP::ClusterMode mode );

  OrbitList(){};
  bool initialize( ODTSP::Problem* prob, ODTSP::ClusterMode mode );

  ~OrbitList(){};
  
  // set functions
  // option a: specify sequence, entry angles, and orbit directions
  // it is assumed that entryAngles and orbitDirs are sorted corresponding
  // to the seq. e.g., entryAngles[i] is assigned to target indicated by seq[i]
  bool set_orbitProperties( vi & seq, vd & entryAngles, 
                            vi & orbitDirs );

  // option b: specify entry states
  bool set_orbitProperties( vi & seq, vd & entryX, 
                            vd & entryY,  vd & entryHrad );

  // option c: specify exit states
  bool set_orbitPropertiesExit( vi & seq, vd & exitX, 
                                vd & exitY, vd & exitHrad );

  // when using maximal clusters tuning
  bool set_tuningParam( double tuningParam ){ 
    tuningParam_ = tuningParam; 
    return true;
  ;}
  
  // get the size of the orbit list 
  int size(){ return list_.size(); };


  // display functions
  void print(); // print all the targets, including any clusters
  void printClusters(); // only print clusters
  
  // get functions
  double get_centerX(int ind);
  double get_centerY(int ind);
  double get_orbitRadius(int ind);
  double get_orbitAngleRad(int ind);
  int get_orbitDir(int ind); 

  double get_entryAngleRad(int ind);
  double get_entryX(int ind);
  double get_entryY(int ind);
  double get_entryHrad(int ind);

  double get_exitAngleRad(int ind);
  double get_exitX(int ind);
  double get_exitY(int ind);
  double get_exitHrad(int ind);


  ODTSP::TargetType get_type(int ind);
  int get_numElems(int ind);
  vi get_elems(int ind);
  bool contains_elem( int listInd, int elemInd);
  bool contains_elem( int listInd, int elemInd, int & entryInd);
  bool elemsInCommon(int i, int j, vi & elems);
  bool elemsInCommon(int i, int j);

  ODTSP::ClusterMode get_mode(){return mode_;};
  std::vector<ODTSP::Target> get_targetList(){return list_;};

}; // class

// Given a n-sized vector 'baseCombo', find the k-sized subsets that contains 
// all combinations of the baseCombo vector.
// e.g., baseCombo: {0,1,2,3}, then n = 4 (size of baseCombo), 
// k = 2 (size of subset)
// Combinations = 0,1 : 0,2 : 0,3 : 1,2 : 1,3 : 2,3
// Inputs: baseCombo, k
// Outputs: comboSet is modified via 'push_back' and the enumerated combos are 
// added 
void enumerateCombinations(const vi & baseCombo, int k, 
                           vvi & comboSet);

// debugging   
void printClusterVector( const vi & cluster );
void printComboSet( const vvi & comboSet );
void printBinnedCombos( const vvvi & binnedCombos );

// returns true if sequences are "equal". Vectors that have the same 
// elements are considered repeats. 
// The input vectors are assumed to be sorted.  
bool equal_seq(vi & first, vi & second);

// returns true if the sorted vector first is contained entirely in the sorted
// vector seocnd
bool subset_seq(vi & first, vi & second);

// Assume that the vectors 'first' and 'second' are already sorted in ascending
// integer order.
// If one vector is a subset of the other, then this subset is smaller
// e.g., first = {1,2,3} and second = {1,2,3,4}, then second > first
// If the vectors are of the same size but differ by one digit, then
// the vector with the larger different digit is greater
// e.g., first = {1,2,3} and second = {1,2,4}, then second > first
// If the vectors are of different size and unique, then the ordering
// is based on the first digit to be different
// e.g., first = {2,3} second = {3,5,6}, then second > first
bool compare_seq(vi & first, vi & second);

// returns true if first has more elements than second
bool compare_seq_longest(vi & first, vi & second);

} //namespace

#endif
