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

#include "ODTSP_Dopt.h"
#include "ODTSP_Utils.h"
// stuff
#include <math.h>       /* cos */
#include <stdlib.h>     /* system */
// reading input .tour file
#include <iostream>
#include <fstream>
#include <string>
#include<limits>
#include<algorithm>
// Artur's Work
#include "RobustDubins_Problem.h"
#include "RobustDubins_Solver.h"
#include "RobustDubins_Path.h"

#include "ConvectedDubins_Problem.h"
//#include<ConvectedDubins_Solver.h>
#include "ConvectedDubins_Path.h"

#include "MathTools.h"

#include <typeinfo> // debugging


ODTSP::Dopt::Dopt(ODTSP::Problem* prob, ODTSP::OrbitList* list, 
                  int numPointsPerCluster, int precision, std::string filename){

  // initialize
  prob_ = prob;
  list_ = list;
  precision_ = precision;
  name_ = filename; 
  M_ =  100;
  offset_ = 0;
  maxVal_ = std::numeric_limits<int>::max()/10;

  convertedFlag = false;  
  ne_ = numPointsPerCluster;
  nt_ = 0;
  nc_ = 0;
  nl_ = 0;
  nr_ = 0;
  nn_ = 0;

}

ODTSP::Dopt::Dopt( ODTSP::DoptSolvers sol, 
                   ODTSP::Problem* prob, ODTSP::OrbitList* list, 
                   int precision, std::string filename){
  if ( sol == ETSP_NO_CLUSTERS ){
    // initialize
    prob_ = prob;
    list_ = list;
    precision_ = precision;
    name_ = filename; 
    sol_ = sol;
    M_ =  std::numeric_limits<int>::max()/2;
    nl_ = list_->size();  
    nn_ = nl_;
    C_ = new int*[nn_];
    for(int i = 0; i < nn_; i++){
      C_[i] = new int[nn_];
        for (int j = 0; j < nn_; j++){
          C_[i][j]=  -precision_;
        }
    }
    buildCostMatrix();
    
    convertedFlag = true;
    writeAsymTspParFile();
    solveATSP();
  }
  else {
    printf("ODTSP::Dopt: Error, incorrect solver specified.\n");
  }
}


ODTSP::Dopt::~Dopt(){
  // clean-up
  for(int i = 0; i < nn_; i++){
      delete [] C_[i];
  }
  delete [] C_;
}

void ODTSP::Dopt::solve(ODTSP::DoptSolvers sol){
  sol_ = sol;
  initialize();
  buildCostMatrix();
  writeGTSPfile();
  if (    sol_ == ODTSP::DoptSolvers::GLKH_NO_CLUSTERS
       || sol_ == ODTSP::DoptSolvers::GLKH_INTERSECTING_CLUSTERS ){
    solveGTSP();
  }
  else if (    sol_ == ODTSP::DoptSolvers::LKH_INTERSECTING_CLUSTERS 
            || sol_ == ODTSP::DoptSolvers::LKH_NO_CLUSTERS ){
    solveATSP();
  }
  else {
    throw std::runtime_error("ODTSP::Dopt::solve: invalid solver specified.");
  }
  extractSoln(); 
  printSoln();
}


void ODTSP::Dopt::resizeSolnVectors( const int & n){
  fromX_.resize(n);
  fromY_.resize(n);
  fromHrad_.resize(n);
  fromCluster_.resize(n);
  fromNode_.resize(n);
  toX_.resize(n);
  toY_.resize(n);
  toHrad_.resize(n);
  toCluster_.resize(n);
  toNode_.resize(n);
}

void ODTSP::Dopt::initialize(){
  nl_ = list_->size();  
  entryAngleRad_.resize(ne_);
  orbitDir_.resize(ne_);  

  // compute total number of nodes, targets, clusters 
  if ( sol_ == GLKH_NO_CLUSTERS || sol_ == LKH_NO_CLUSTERS ){
    nn_ = nl_*ne_;
    nt_ = nl_;  
    resizeSolnVectors(nl_);
  }
  else if (    sol_ == GLKH_INTERSECTING_CLUSTERS 
            || sol_ == LKH_INTERSECTING_CLUSTERS ){
    // go through entire list 
    for ( int i = 0; i < nl_; i++){
      if ( list_->get_type(i)==ODTSP::TargetType::Cluster ){
        nc_++; // increment number of clusters 
      }
      else if ( list_->get_type(i)==ODTSP::TargetType::Single ){
        nt_++; // increment number of targets 
      }
    }
    // go through the clusters 
    for (int i = 0; i < nc_; i++){
        // get vector of indices of the targets in the i-th cluster           
        vi elems = list_->get_elems(nt_ + i);
        // also, append them to one long vector  
        targetsInClusters_.insert( targetsInClusters_.end(), 
                                   elems.begin(),
                                   elems.end() );
        nr_ += list_->get_numElems(nt_ + i) - 1; // increment the number of replicas required 
    }
    nn_ = ne_*(nt_ + nc_+ nr_); // total number of nodes
     // targetsInClusters_ can have repetitive entries so sort and remove duplicates
    std::sort(targetsInClusters_.begin(), targetsInClusters_.end());
    vi::iterator it;
    it = std::unique(targetsInClusters_.begin(), targetsInClusters_.end() );
    targetsInClusters_.resize( std::distance(targetsInClusters_.begin(),it) );
  
    printGraphProperties();
    // popule node list index which assigns to each node the corresponding
    // entry in the list_
    nodeListIndex_.resize(nn_);
    clusterSize_.resize(nc_);
    pastReplicas_.resize(nc_);
    for (int i = 0; i < nt_; i++){
      for (int j = 0; j < ne_; j++){
        nodeListIndex_[i*ne_+j]=i;
      }
    }
    int qxSum = 0;
    for (int i = 0; i < nc_; i++){
      pastReplicas_[i] = qxSum;
      clusterSize_[i] = list_->get_numElems(nt_+i);
      for (int j = 0; j < clusterSize_[i]; j++){
        for (int k = 0; k < ne_; k++){
          nodeListIndex_[ne_*(nt_+qxSum)+k]=nt_+i;
        }
        qxSum++;
      }
    }
    //printf("pastReplicas_\n");
    //printClusterVector(pastReplicas_);
  }
  // allocate memory for the cost matrix C_
  // initialize all values of C_ to -precision_
  // (this value indicates undefined and causes printing `?' entry to screen)
  C_ = new int*[nn_];
  for(int i = 0; i < nn_; i++){
    C_[i] = new int[nn_];
      for (int j = 0; j < nn_; j++){
        C_[i][j]=  -precision_;
      }
  }
  // generate the entry angles starting from 0 (along + x-axis, going CCW)
  for(int i = 0; i < ne_ / 2; ++i){
    entryAngleRad_[i*2]   = (double) ( ((double)i)/((double)ne_/2.0) * 2.0*M_PI );
    entryAngleRad_[i*2+1] = entryAngleRad_[i*2];
    orbitDir_[i*2]   = -1; // CW orbit
    orbitDir_[i*2+1] = 1; // CCW orbit
  }
}

void ODTSP::Dopt::buildCostMatrix(){

  if ( sol_ == ETSP_NO_CLUSTERS ){
    // populate the initial cost matrix, here we ignore the cycles for a minute. 
    for(int i = 0; i < nl_; i++){ // from orbit i
      for(int j = 0; j < nl_; j++){ // to orbit j 
        double dx = list_->get_centerX(j) -  list_->get_centerX(i);
        double dy = list_->get_centerY(j) -  list_->get_centerY(i);
        C_[i][j] = (int)((dx*dx + dy*dy)*(double)precision_);
      }
    }
  }
  else if ( sol_ == GLKH_NO_CLUSTERS ){
    // populate the initial cost matrix, here we ignore the cycles for a minute. 
    for(int i = 0; i < nl_; i++){ // from orbit i
      for(int j = 0; j < nl_; j++){ // to orbit j
        for(int k = 0; k < ne_; k++){ // from k-th node of orbit i
	        for(int l = 0; l < ne_; l++){ // to l-th node of orbit j 
            // disconnect nodes from themself by giving 'infinite' cost
	          if(i == j){
	            C_[ne_*i+k][ne_*j+l] = maxVal_;//std::numeric_limits<int>::max();
            }	        
            else { // compute the integer Dubins path 
	            C_[ne_*i+k][ne_*j+l] = dubinsCost(i,k,j,l);
	          }
	        }
        }
      }
    }
  }
  else if ( sol_ == LKH_NO_CLUSTERS ){
    // populate the initial cost matrix, here we ignore the cycles for a minute. 
    for(int i = 0; i < nl_; i++){ // from orbit i
      for(int j = 0; j < nl_; j++){ // to orbit j
        // note: ne_ is twice the number of entry angles (includes direction)
        for(int k = 0; k < ne_; k++){ // from k-th node of orbit i
	        for(int l = 0; l < ne_; l++){ // to l-th node of orbit j 
            // disconnect nodes from themself by giving 'infinite' cost
	          if(i == j){
              if ( k == l-1 ){ // create a zero-cost cycle 
                C_[ne_*i+k][ne_*j+l] = 0;
              }
              else if ( k == ne_-1 && l == 0){ 
                C_[ne_*i+k][ne_*j+l] = 0;
              }
              else { 
                C_[ne_*i+k][ne_*j+l] = maxVal_; 
              }
            }	        
            else { // compute the integer Dubins path 
              if ( k == ne_ ){
                C_[ne_*i+k][ne_*j+l] = dubinsCost(i,0,j,l) + M_;
              }
              else{
	              C_[ne_*i+k][ne_*j+l] = dubinsCost(i,k+1,j,l) + M_;
              }
	          }
	        }
        }
      }
    }
  }
  else if ( sol_ == GLKH_INTERSECTING_CLUSTERS ){
    offset_ = 0;
    // populate the initial cost matrix, here we ignore the cycles for a minute. 
    for(int i = 0; i < nt_; i++){ // from orbit i
      for(int j = 0; j < nt_; j++){ // to orbit j
        for(int k = 0; k < ne_; k++){ // from k-th node of orbit i
	        for(int l = 0; l < ne_; l++){ // to l-th node of orbit j 
            // disconnect nodes from themself by giving 'infinite' (max) cost
	          if( i==j ){
	            C_[ne_*i+k][ne_*j+l] = maxVal_;
            }	        
            else { // compute the integer Dubins path cost 
              C_[ne_*i+k][ne_*j+l] = dubinsCost(i,k,j,l) + offset_;
	          }      
	        }
        }
      }
    }
    //printCostMatrix(); 
    //printf("clusterSize_");
    //printClusterVector(clusterSize_);
    int qxCurSum = 0; // total number of replica targets added so far 
    for (int k = 0; k < nc_; k++){ // for all of the clusters 
      for (int i = 0; i < clusterSize_[k]; i++){ // from each target in current cluster 
        for (int j = 0; j < nt_ ; j++){ // to all targets  
          for (int l = 0; l < ne_; l++){ // replica entry angle no. 
            for (int m = 0; m < ne_; m++){ // destination entry angle no. 
              int fromNode = ne_*(nt_ + qxCurSum) + l;
              int toNode = ne_*j + m;
              // disconnect nodes from themself by giving 'infinite' (max) cost
	            if ( list_->contains_elem(nt_+k, j) ){
	                C_[fromNode][toNode] = maxVal_; 
                  C_[toNode][fromNode] = maxVal_; 
              }     
              else { // compute the integer Dubins path cost 
                if ( i == 0 ){ // from first node in cluster 
	                C_[fromNode][toNode] = maxVal_;
                  C_[toNode][fromNode] = dubinsCost(j,m,nt_+k,l) + offset_;
                }
                else if ( i != clusterSize_[k]-1 ){ // middle node in cluster
	                C_[fromNode][toNode] = maxVal_; 
                  C_[toNode][fromNode] = maxVal_;
                }
                else { // final node 
	                C_[fromNode][toNode] = dubinsCost(nt_+k,l,j,m) + offset_;
                  C_[toNode][fromNode] = maxVal_;
                }
	            }
            }
          }
        }
      qxCurSum++;
      }
    }
    //printCostMatrix(); 

    // construct inter-cluster edges 
    int prevReplicas = 0; // replicas added from previous clusters 
    for (int i = 0; i < nc_; i++){ 
      // using j,k,l,m we go through all combinations of repl. nodes in cluster
      for (int j = 0; j < clusterSize_[i]; j++){
        for (int k = 0; k < clusterSize_[i]; k++){
          for (int l = 0; l < ne_; l++){   
            for (int m = 0; m < ne_; m++){   
              int fromNode = ne_*(nt_ + prevReplicas + j) + l;
              int toNode = ne_*(nt_ + prevReplicas + k) + m;  
              if ( k == j+1 && l == m ){  
                // create zero-cost directed edge to succesive node 
                // representing same entry angle 
                  C_[fromNode][toNode] = offset_;             
              } // otherwise disconnect inter-cluster nodes 
              else { 
                C_[fromNode][toNode] = maxVal_;
              }
            }
          }
        }       
      }  
      prevReplicas += clusterSize_[i];   
    }

    // construct intra-cluster edges 
    int qxk_sum = 0; // total number of replica targets added so far 
    for (int k = 0; k < nc_; k++){ // number of clusters 
      int qxp_sum = 0;
      for ( int p = 0; p < nc_; p++){
        if ( p != k ){
          // consider all combinations of i and j targets in k-th cluster
          for (int i = 0; i < clusterSize_[k]; i++){
            for ( int j = 0; j < clusterSize_[p]; j++){
              // for all combinations of entry angles 
              for (int l = 0; l < ne_; l++){   
                for (int m = 0; m < ne_; m++){   
                  int fromNode = ne_*(nt_ + qxk_sum + i) + l;
                  int toNode = ne_*(nt_ + qxp_sum + j) + m;
                  // check if clusters have no common elements
                  if ( !(list_->elemsInCommon(nt_+k,nt_+p) )
                       && i == clusterSize_[k]-1 && j == 0){  
                    // (i == clusterSize_[k]-1) means only connect the last 
                    // set of replica nodes 
                    C_[fromNode][toNode] = dubinsCost(nt_+k,l,nt_+p,m) + offset_;          
                  }
                  else {
                    C_[fromNode][toNode] = maxVal_;
                  }
                }
              }
            }     
          }          
        }  
        qxp_sum += clusterSize_[p];
      }
      qxk_sum += clusterSize_[k];  
    }

  }
  else if ( sol_ == LKH_INTERSECTING_CLUSTERS ){
    // compute Sx_ 
    for ( int i = 0; i < list_->size(); i++ ){
      Sx_.push_back( list_->get_elems(i) );
    }

    // populate the initial cost matrix, here we ignore the cycles for a minute. 
    for(int i = 0; i < nt_; i++){ // from orbit i
      for(int j = 0; j < nt_; j++){ // to orbit j
        // check if one of the orbits is a cluster
        bool selfConnected = false;
        if ( i==j ){
          selfConnected = true;
        }        
        for(int k = 0; k < ne_; k++){ // from k-th node of orbit i
	        for(int l = 0; l < ne_; l++){ // to l-th node of orbit j 
            // disconnect nodes from themself by giving 'infinite' (max) cost
	          if( selfConnected ){
	            C_[ne_*i+k][ne_*j+l] = std::numeric_limits<int>::max();
            }	        
            else { // compute the integer Dubins path cost 
              if ( k == ne_ && l == ne_ ){
                C_[ne_*i+k][ne_*j+l] = dubinsCost(i,0,j,0);
              }
              if ( k == ne_ ){
                C_[ne_*i+k][ne_*j+l] = dubinsCost(i,0,j,l+1);
              }
              if ( l == ne_ ){
                C_[ne_*i+k][ne_*j+l] = dubinsCost(i,k+1,j,0);
              }
              else{
	              C_[ne_*i+k][ne_*j+l] = dubinsCost(i,k+1,j,l+1);
              }
	          }      
	        }
        }
      }
    }
    // go through and create the cycles in each group. 
    // There are nl groups with ne points. There should a cyle ne long.
    for(int i = 0; i < nt_; i++){
      for(int j = 0; j < ne_; j++) {
        // if we're at the last vertex in the tour, and it is not in a cluster  
        if( j == ne_ - 1 && std::find( targetsInClusters_.begin(),
                                       targetsInClusters_.end(), i ) == 
                                       targetsInClusters_.end() ) { 
	        C_[ne_*i+j][ne_*i] = 0;
        } 
        if ( j != ne_ - 1) {
	        C_[ne_*i+j][ne_*i+j+1] = 0;
        }	
      }
    }
    // add M_ to all inter-cluster edges (do not add to intra-cluster edges)
    for(int i = 0; i < nt_; i++){
      for(int j = 0; j < nt_; j++){
        for(int k = 0; k < ne_; k++){
	        for(int l = 0; l < ne_; l++){
            // if already disconnected to do not add M_ this will cause overflow   
            if(  i != j && 
                      C_[ne_*i+k][ne_*j+l] != std::numeric_limits<int>::max() ){
	            C_[ne_*i+k][ne_*j+l] += M_;
            }
	        }
        }
      }
    }
    // we now create zero-cost edges connecting replica nodes in Sx_:
    //    - the first replicas begin after ne_*nt_ rows 
    //    - if the k-th target appears in qx clusters then there will be 
    //      a total of ne_*qx replicas for that Sx_[k]
    //    - by definition of Sx_, qx >= 2
    //    - replicas of, in general, begins after row ne_*nt_ + ne_*qxCurSum
    // for each target 
    int qxCurSum = 0; // total number of replica targets added so far 
    for (int k = 0; k < nc_; k++){ // number of clusters 
      int qx = Sx_[k].size(); // number of targets in the k-th cluster             
      // consider all combinations of i and j targets in k-th cluster
      for (int i = 0; i < qx; i++){
        for ( int j = 0; j < qx; j++){
          // for all combinations of entry angles 
          for (int l = 0; l < ne_; l++){   
            for (int m = 0; m < ne_; m++){   
              int fromNode = ne_*(nt_ + qxCurSum + i) + l;
              int toNode = ne_*(nt_ + qxCurSum + j) + m;  
              if ( i == j ){ // within a replica createa cycle 
                if ( l==m-1 ){
                  C_[fromNode][toNode] = 0;
                }
                else {
                  C_[fromNode][toNode] = std::numeric_limits<int>::max();
                }
              } 
              else { // create a zero cost edge between replicas ux,i and ux,j
                if ( (l==m) ){
                  C_[fromNode][toNode] = 0;
                }
                else {
                  C_[fromNode][toNode] = std::numeric_limits<int>::max();
                }
              }
            }
          }
        }
      }   
      qxCurSum++; // a set of ne_ replicas processed     
    }
    qxCurSum = 0; // total number of replica targets added so far 
    // Sx_ set of all indices containing cluster x
    // Sx_[i] gives indices that targetsInClusters_[i] is include

    // Sx_.size() is the number of targets part of clusters 
    for (int k = 0; k < Sx_.size(); k++){ // for each target in a cluster
      vi clusterElems = Sx_[k];
      int qx = Sx_[k].size(); // no. of clusters the k-th target is part of 
      for (int i = 0; i < qx; i++){ // for each cluster in  replica no. in qx 
        for (int j = 0; j < nt_ ; j++){ // to target/replica no.s
          for (int l = 0; l < ne_; l++){ // replica entry angle no. 
            for (int m = 0; m < ne_; m++){ // destination entry angle no. 
              bool selfConnected = false;
//              if ( std::find( clusterElems.begin(),
//                              clusterElems.end(), j ) != 
//                              clusterElems.end() ){
//                selfConnected = true;                     
//              }
              if ( clusterElems[i]==j ){
                selfConnected = true;
              }
              int fromNode = ne_*(nt_ + qxCurSum + i) + l;
              int toNode = ne_*j + m;
              // disconnect nodes from themself by giving 'infinite' (max) cost
	            if( selfConnected ){
                if ( clusterElems[i]==j && l == ne_-1 && m == 0){ // forms a cycle 
                  C_[fromNode][toNode] = 0;
                  C_[toNode][fromNode] = std::numeric_limits<int>::max();  
                }
                else if ( selfConnected && l == 0 && m ==ne_-1 ){
                  C_[fromNode][toNode] = std::numeric_limits<int>::max();
                  C_[toNode][fromNode] = 0;  
                }
                else {
	                C_[fromNode][toNode] = std::numeric_limits<int>::max(); 
                  C_[toNode][fromNode] = std::numeric_limits<int>::max(); 
                }

              }	        
              else { // compute the integer Dubins path cost 
	              C_[fromNode][toNode] = dubinsCost(nt_+k,l,j,m) + M_;
                C_[toNode][fromNode] = dubinsCost(j,m,nt_+k,l) + M_;
	            }
            }
          }
        }
      }
    }
  }
}


void ODTSP::Dopt::plotCostMatrixMatlab(){
  std::string filename = "tempMatrix.m";
  FILE * pFile = fopen(filename.c_str(),"w");
  for (int i=0; i < nn_; i++){
    for (int j=0; j < nn_; j++){
      fprintf(pFile,"C(%i,%i)=%i;",i+1,j+1,C_[i][j]);
    }
  }
  fclose(pFile);
}

int ODTSP::Dopt::dubinsCost( const int & fromTargInd,
                             const int & fromNodeInd, 
                             const int & toTargInd,
                             const int & toNodeInd ){
  // temp variables
  double fromX, fromY, fromHrad;
  double toX, toY, toHrad;
  double exitAngle, dir;

  // get the "from" point (i-th cluster, k-th node) + orbitAngle
  // (the x,y coordiantes are NOT scaled up)
  dir = (double)orbitDir_[fromNodeInd];
  exitAngle = entryAngleRad_[fromNodeInd] 
              + dir*(list_->get_orbitAngleRad(fromTargInd)); 

  ODTSP::transformEntryAngle(list_->get_orbitRadius(fromTargInd), 
                             exitAngle, dir, 
                             list_->get_centerX(fromTargInd), 
                             list_->get_centerY(fromTargInd), 
                             // output
                             fromX, fromY, fromHrad);

  // get the "to" point (j-th cluster, l-th node)
  dir = (double)orbitDir_[toNodeInd];
  ODTSP::transformEntryAngle(list_->get_orbitRadius(toTargInd), 
                             entryAngleRad_[toNodeInd], dir, 
                             list_->get_centerX(toTargInd), 
                             list_->get_centerY(toTargInd), 
                             // output
                             toX, toY, toHrad);
  
  int cost;  
  // define the Dubins problem and solve
  if ( prob_->motionModel() == DUBINS){
    RobustDubins::Problem problemStatement;	        
    problemStatement.set_minTurningRadius(prob_->R() );
    problemStatement.set_stateInitial(fromX, fromY, fromHrad);
    problemStatement.set_stateFinal(toX, toY, toHrad);

    RobustDubins::Solver rds;  
    rds.set_problemStatement(problemStatement);
    rds.solve();
    cost = (int) (rds.get_optCost() * precision_);
  }
  else if ( prob_->motionModel() == CONVECTED_DUBINS ){
    ConvectedDubins::Problem cDubinsProb;
    cDubinsProb.set_vehicleProperties( prob_->R() , prob_->v() );
    // TODO: Change wind direction, 5-Jan-2018         
    //cDubinsProb.set_windMagDir( prob_->w() , prob_->hwRad() );
    cDubinsProb.set_windMagDir( prob_->w() , prob_->hwRadNED() );
    cDubinsProb.set_stateInitial_XYEastNorthCourse( fromX, 
                                                    fromY, 
                                                    fromHrad );
    cDubinsProb.set_stateFinal_XYEastNorthCourse( toX, 
                                                  toY, 
                                                  toHrad );
    //cDubinsProb.print();
    ConvectedDubins::Solver cds( &cDubinsProb );
    cost = (int) (cds.get_optCost() * precision_);    
  }
  else if ( prob_->motionModel() == FEASIBLE_DUBINS ){
    RobustDubins::Problem problemStatement;
    double eps = prob_->w() / prob_->v();
    double feasibleR = prob_->R()*(1.0 + eps)*(1.0 + eps);         
    problemStatement.set_minTurningRadius( feasibleR );
    problemStatement.set_stateInitial(fromX, fromY, fromHrad);
    problemStatement.set_stateFinal(toX, toY, toHrad);

    RobustDubins::Solver rds;  
    rds.set_problemStatement(problemStatement);
    rds.solve();
    // input for optPathTimeInCurrents requires current in ENU not NED 
    cost = (int) (rds.get_optPathTimeInCurrents( prob_->v(), prob_->w(), prob_->hwRad() ) * precision_);
  }
  if ( cost < 0 ){
    throw std::runtime_error("ODTSP::Dopt::dubinsCost < 0\n"); 
  }
  cost = cost + orbitCost( toTargInd , toNodeInd );
  return cost;
}

int ODTSP::Dopt::orbitCost(  const int & toTargInd,
                             const int & toNodeInd ){

  // both dubins and convected dubins need this data:
  double orbitRadius = list_->get_orbitRadius(toTargInd);
  double orbitAngleRad = list_->get_orbitAngleRad(toTargInd);
  int cost;
   
  // define the Dubins problem and solve
  if ( prob_->motionModel() == DUBINS){
    // then orbit cost is length 
    double orbitLength = orbitRadius*orbitAngleRad;
    cost = (int) ( orbitLength * precision_);
  }
  else if ( prob_->motionModel() == CONVECTED_DUBINS || prob_->motionModel() == FEASIBLE_DUBINS ){
    // get the "from" point (i-th cluster, k-th node) + orbitAngle
    // (the x,y coordiantes are NOT scaled up)
    double dir = (double)orbitDir_[toNodeInd];
    double toX, toY, toHrad;
    ODTSP::transformEntryAngle(list_->get_orbitRadius(toTargInd), 
                               entryAngleRad_[toNodeInd], dir, 
                               list_->get_centerX(toTargInd), 
                               list_->get_centerY(toTargInd), 
                               // output
                               toX, toY, toHrad);
    // not that the above transformEntryAngle/orbitlist function uses ENU frame
    // convection of X - East , angles +CCW from East 
    // dir +1 implies CCW 

    // but the below computation derived in paper assumes NED frame 
    // convention of X - North, angles +CW from North
    // dir +1 implies CW 

    // convert to XY = NE coordinate system 
    double v = prob_->v();
    double eps = prob_->w()/v; 
    double chiNED = fmod( M_PI/2.0 - toHrad , 2.0*M_PI );
    double varphiNED = prob_->hwRadNED(); // psi is measured from X axis north CW
    double psiNED;

    
    // convert dir
    double dirMod = -dir;    
    int N = 100;
    double delth = prob_->get_orbitAngleRad()/( double(N) );
    double timeSum = 0;
    
      
    double vg, delt;
    for ( int i = 0; i < N; i ++ ){

      psiNED = fmod( chiNED - asin( eps*sin(varphiNED - chiNED) ), 2.0*M_PI );
      vg = v*sqrt(1.0 + eps*( eps + 2.0*cos(varphiNED - psiNED) ) ) ;   
      delt = orbitRadius/vg*delth;
      // update 
      timeSum = timeSum + delt;
      chiNED = chiNED + delth*dirMod;
    }
    cost = (int) ( timeSum * precision_);
  }
  if ( cost < 0 ){
    throw std::runtime_error("ODTSP::Dopt::orbitCost < 0\n"); 
  }
  return cost;
}


void ODTSP::Dopt::writeGTSPfile(){
  // set up the file name
  std::string fName = name_;
  fName.append(".gtsp");
  printGraphProperties();
  if ( sol_ == LKH_INTERSECTING_CLUSTERS || sol_ == LKH_NO_CLUSTERS){
    FILE * pFile;
    pFile = fopen (fName.c_str(),"w");
    fprintf(pFile,"NAME:  dubins\n");
    fprintf(pFile,"TYPE: ATSP\n");
    fprintf(pFile,"COMMENT: %d city problem\n", nn_);
    fprintf(pFile,"DIMENSION:  %d\n", nn_);
    fprintf(pFile,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(pFile,"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    // print the full matrix to define the gtsp 
    fprintf(pFile,"EDGE_WEIGHT_SECTION\n");
    for(int i = 0; i < nn_; i++){
      for(int j = 0; j < nn_; j++){
	      fprintf(pFile," %d ", C_[i][j]);
      }
      fprintf(pFile,"\n");    
    }
    fprintf(pFile,"EOF\n");
    fclose(pFile);
  }
  else if ( sol_ == GLKH_INTERSECTING_CLUSTERS ){
    FILE * pFile;
    pFile = fopen (fName.c_str(),"w");
    fprintf(pFile,"NAME:  dubins\n");
    fprintf(pFile,"TYPE: AGTSP\n");
    fprintf(pFile,"COMMENT: %d city problem\n", nn_);
    fprintf(pFile,"DIMENSION:  %d\n", nn_);
    // New entry GTSP_SETS: M (where M is a positive integer) defines the 
    // number of clusters in this instance. This is used in conjunction with 
    // GTSP_SET_SECTION defined below. 
    fprintf(pFile,"GTSP_SETS: %d\n", nt_);
    fprintf(pFile,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(pFile,"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    // print the full matrix to define the gtsp 
    fprintf(pFile,"EDGE_WEIGHT_SECTION\n");
    for(int i = 0; i < nn_; ++i){
      for(int j = 0; j < nn_; ++j){
	      fprintf(pFile," %d ", C_[i][j]);
      }
      fprintf(pFile,"\n");    
    }
    // GTSP_SET_SECTION defines which nodes belong to which clusters
    // This section contains exactly M entries, where M is the number of clusters 
    // Each entry has the following format:
    // m v1 v2 ... vk(m) -1, where m is the cluster number 
    // (clusters are numbered from 1 to M), and 
    // v1 v2 ... vk(m) are vertices comprising cluster m (vertices from 1 to N). 
    // http://www.cs.rhul.ac.uk/home/zvero/GTSPLIB/
    fprintf(pFile,"GTSP_SET_SECTION:\n");
    for(int i = 0; i < nt_; i++) {
      //printf(" ***ADDING NEW SECTION FOR TARGET %i \n",i);
      fprintf(pFile,"%d",i+1);
      // add the nodes forming an orbit around i-th target:
      for(int j = 0; j < ne_; j++){
        fprintf(pFile," %d",i*ne_+j+1);
      }
      // add additional replica nodes (if any)
      for (int j = 0; j < nc_; j++){
        int entryInd; 
        //printf(" j = %i, pastReplicas = %i \n",j, pastReplicas_[j]);
        if ( list_->contains_elem(nt_+j, i, entryInd) ){
          //printf(" cluster no. %i contains target %i in the %i-th position \n",j,i,entryInd);
          //printf(" list_->contains_elem(nt_ + j, i, entryInd) = (%i,%i,%i)\n",nt_+j,i,entryInd);
          // the i-th target is in cluster j, then add the ne_ replica nodes 
          // corresponding to the position in the cluster 
          for (int k = 0; k < ne_; k++){
            // the 1 is added because GTSP format starts with node = 1 not 0
            fprintf(pFile," %d", ne_*(nt_+pastReplicas_[j]+entryInd)+k+1);
            //printf(" %d", ne_*(nt_+pastReplicas_[j]+entryInd)+k+1);
          }
          //printf("\n");
        }
      }
      fprintf(pFile," -1\n");
      //printf("-1 \n");
    }
    fprintf(pFile,"EOF\n");
    fclose(pFile);
  }
  else if ( sol_ == GLKH_NO_CLUSTERS ){
    FILE * pFile;
    pFile = fopen (fName.c_str(),"w");
    fprintf(pFile,"NAME:  dubins\n");
    fprintf(pFile,"TYPE: AGTSP\n");
    fprintf(pFile,"COMMENT: %d city problem\n",nl_ * ne_);
    fprintf(pFile,"DIMENSION:  %d\n",nl_ * ne_);
    fprintf(pFile,"GTSP_SETS: %d\n",nl_);
    fprintf(pFile,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(pFile,"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    // print the full matrix to define the gtsp 
    fprintf(pFile,"EDGE_WEIGHT_SECTION\n");
    for(int i = 0; i < nl_*ne_; ++i){
      for(int j = 0; j < nl_*ne_; ++j){
	      fprintf(pFile," %d ", C_[i][j]);
      }
      fprintf(pFile,"\n");    
    }
    // GTSP_SET_SECTION defines which nodes belong to which clusters
    // This section contains exactly M entries, where M is the number of clusters 
    // Each entry has the following format:
    // m v1 v2 ... vk(m) -1, where m is the cluster number 
    // (clusters are numbered from 1 to M), and 
    // v1 v2 ... vk(m) are vertices comprising cluster m (vertices from 1 to N). 
    // http://www.cs.rhul.ac.uk/home/zvero/GTSPLIB/
    fprintf(pFile,"GTSP_SET_SECTION:\n");
    for(int i = 1; i < nl_+1; ++i) {
      fprintf(pFile,"%d",i);
      for(int j = 1; j < ne_+1; ++j){
        fprintf(pFile," %d",(i-1)*ne_+j);
      }
      fprintf(pFile," -1\n");
    }
    fprintf(pFile,"EOF\n");
    fclose(pFile);
  }


  // write par file 
  // For more info see LKH-2_USERS_GUIDE.pdf (July 09)
  FILE * pFile2;
  std::string fName2 = name_;
  fName2.append(".par");
  pFile2 = fopen (fName2.c_str(),"w");
  fprintf(pFile2,"PROBLEM_FILE = %s.gtsp\n",name_.c_str());    

  //fprintf(pFile2,"ASCENT_CANDIDATES = 500\n");
  //fprintf(pFile2,"INITIAL_PERIOD = 1000\n");

  //fprintf(pFile2,"MAX_TRIALS = 1000\n");
  // fprintf(pFile2,"MAX_CANDIDATES = 5\n");
  // MAX_SWAMPS
  // MERGE_TOUR_FILE
  // MOVE_TYPE
  // NONSEQUENTIAL_MOVE_TYPE
  fprintf(pFile2,"OUTPUT_TOUR_FILE = %s.tour\n",name_.c_str());
  // OPTIMUM 
  // PATCHING_A
  // PATCHING_C
  fprintf(pFile2,"PI_FILE = %s.pi\n",name_.c_str());
  //fprintf(pFile2,"POPULATION_SIZE = 0\n");
  fprintf(pFile2,"PRECISION = 1\n");
  // RESTRICTED SEARCH
  //fprintf(pFile2,"RUNS = 10\n");
  //fprintf(pFile2,"SEED = 1\n");
  // STOP_AT_OPTIMUM
  // SUBGRADIENT
  // SUBPROLBME_SIZE
  // SUBPROBLEM_TOUR_FLIGHT 
  // SUBSEQUENT_MOVE_TYPE
  // SUBSEQUENT_PATCHING 
  // TIME_LIMIT 
  fprintf(pFile2,"TOUR_FILE = %s.tour\n",name_.c_str());  
  fprintf(pFile2,"TRACE_LEVEL = 1\n"); // has to do with verbosity
  fclose(pFile2);
}


void ODTSP::Dopt::writeAsymTspParFile(){
  // set up the file name
  std::string fName = name_;
  fName.append(".tsp");

  // write the .tsp and .par files
  if(convertedFlag == false)
    printf("ATSP file not printed, not converted to ATSP!\n");
  else {
    // write the tsp file
    FILE * pFile;
    pFile = fopen (fName.c_str(),"w");
    fprintf(pFile,"NAME:  dubins\n");
    fprintf(pFile,"TYPE: ATSP\n");
    fprintf(pFile,"COMMENT: %d city problem\n",nl_ * ne_);
    fprintf(pFile,"DIMENSION:  %d\n",nn_);
    fprintf(pFile,"EDGE_WEIGHT_TYPE: EXPLICIT\n");
    fprintf(pFile,"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n");
    fprintf(pFile,"EDGE_WEIGHT_SECTION\n");
    for(int i = 0; i < nn_; ++i){
      for(int j = 0; j < nn_; ++j){      
	      fprintf(pFile," %d ", C_[i][j]);
      }
      fprintf(pFile,"\n");    
    }
    fprintf(pFile,"EOF\n");
    fclose(pFile);

    // write par file 
    FILE * pFile2;
    std::string fName2 = name_;
    fName2.append(".par");
    pFile2 = fopen (fName2.c_str(),"w");
    fprintf(pFile2,"PROBLEM_FILE = %s.tsp\n",name_.c_str());    
    fprintf(pFile2,"TRACE_LEVEL = 0\n");
    fprintf(pFile2,"TOUR_FILE = %s.tour\n",name_.c_str());
    fprintf(pFile2,"ASCENT_CANDIDATES = 500\n");
    fprintf(pFile2,"INITIAL_PERIOD = 1000\n");
    fprintf(pFile2,"MAX_CANDIDATES = 30\n");
    fprintf(pFile2,"MAX_TRIALS = 1000\n");
    fprintf(pFile2,"OUTPUT_TOUR_FILE = %s.tour\n",name_.c_str());
    fprintf(pFile2,"PI_FILE = %s.pi\n",name_.c_str());
    fprintf(pFile2,"POPULATION_SIZE = 0\n");
    fprintf(pFile2,"PRECISION = 1\n");
    fprintf(pFile2,"RUNS = 10\n");
    fprintf(pFile2,"SEED = 1\n");
    fprintf(pFile2,"TRACE_LEVEL = 1\n");
    fclose(pFile2);

  }
}


void ODTSP::Dopt::solveGTSP(){
  // solve the GTSP
  if( sol_ == ODTSP::DoptSolvers::GLKH_NO_CLUSTERS || 
      sol_ == ODTSP::DoptSolvers::GLKH_INTERSECTING_CLUSTERS ) { 
    // run GLKH
    // **** GLKH needs a ./TMP folder!!!
    // **** GLKH needs a symlink to ./LKH folder!!!
    std::string command = "GLKH ";  
    command.append(name_);
    command.append(".par");
    system(command.c_str()); 

    // read in the text file
    std::string iFile = name_;
    iFile.append(".tour");
    std::string line;
    std:: ifstream myfile (iFile.c_str());
    if(myfile.is_open()){
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    std::getline (myfile,line);
	    while(line.compare("-1") != 0){
        // GLKH starts counting nodes from 1, we count from 0, so subtract 1
        atspSoln_.push_back(atoi(line.c_str())-1);
        std::getline (myfile,line);
      }
	    myfile.close();
      }
    else {
      printf(" unable to open file\n");
    }
  } 
  else {
    printf(" did not specify a valid solver\n");
  }

}

void ODTSP::Dopt::solveATSP(){  
  // write the command to execute LKH
  std::string command = "LKH ";
  command.append(name_);
  command.append(".par");
  printf("****DEBUG****: Running Command: %s \n",command);
  system(command.c_str());

  
  // read in the resulting text file, saving the cost and the tour (etspSoln)
  std::string iFile = name_;
  iFile.append(".tour");
  std::string line; // temporary hold the current line
  std:: ifstream myfile (iFile.c_str());
  atspSoln_.clear(); // stores the solution
  if(myfile.is_open()) {
    std::getline (myfile,line);
    std::getline (myfile,line);
    //etspCost_ = stoi( line.substr(19,line.length()-1) );
    std::getline (myfile,line);
    std::getline (myfile,line);
    std::getline (myfile,line);
    std::getline (myfile,line);
    std::getline (myfile,line);
    while(line.compare("-1") != 0){
      printf( "%s \n", line.c_str() );
	    atspSoln_.push_back(atoi(line.c_str())-1); // adjust node number 
	    std::getline (myfile,line);
	  }
    myfile.close();
  }
}

void ODTSP::Dopt::extractSoln(){    

  if ( sol_ == LKH_NO_CLUSTERS ){
    // warning! We can't guarantee that the solution starts at the first returned node
    gtspSoln_.clear();
    // create a verbose flag to print this to terminal  
    printf("ATSP Solution: \n");
    for(int i = 0; i < atspSoln_.size() -1 ; i++){
      int cost = C_[atspSoln_[i]][atspSoln_[i+1]];
      if ( cost == maxVal_ ){
        cost = -1;  
      }
      else if (cost >= M_ ) {
        cost = cost - M_;
        gtspSoln_.push_back(atspSoln_[i+1]);
      }
      printf(" %d ---(Cost: %d)--> ",atspSoln_[i], cost );
    }
    int cost = C_[atspSoln_[atspSoln_.size()-1]][atspSoln_[0]];
    if ( cost == maxVal_ ){
      cost = -1;  
    }
    else if (cost >= M_ ) {
      cost = cost - M_;
      gtspSoln_.push_back(atspSoln_[0]);    
    }   
    printf(" %d ---(Cost: %d)--> %d \n", atspSoln_[atspSoln_.size()-1], cost, atspSoln_[0]  );
    printf(" \n GTSP Solution:\n");
    for(int i = 0; i < gtspSoln_.size()-1 ; i++){
      printf(" %d --> ",gtspSoln_[i]);
    }
    printf(" %d --> %d ",gtspSoln_[gtspSoln_.size()-1], gtspSoln_[0]);
  }
  else { 
    // create a verbose flag to print this to terminal  
    printf("Direct Solution: \n");
    for(int i = 0; i < atspSoln_.size() -1 ; i++){
      printf(" %d ---(Cost: %d)--> ",atspSoln_[i], C_[atspSoln_[i]][atspSoln_[i+1]] );
    }
    printf(" %d ---(Cost: %d)--> ",atspSoln_[atspSoln_.size()-1], C_[atspSoln_[atspSoln_.size()-1]][atspSoln_[0]] );
      
    gtspSoln_.clear();    
    // Create a verbose flag to print this to the terminal  
    printf("\n GTSP Solution: \n");
    int pastCluster  = -1;
    for(int i = 0; i < atspSoln_.size(); i++){
      if ( isReplica( atspSoln_[i] ) ){
        int curCluster = whichCluster( atspSoln_[i] );
        printf(" node %i is in cluster %i  and is replica \n", atspSoln_[i], curCluster);
        if ( curCluster != pastCluster ){
          gtspSoln_.push_back( atspSoln_[i] );
          printf(" %dR --> ", gtspSoln_[ gtspSoln_.size()-1 ] );
          pastCluster = curCluster;
        }
      }
      else {
        printf(" node %i is not a replica \n", atspSoln_[i]);
        gtspSoln_.push_back( atspSoln_[i] );
        printf(" %d --> ", gtspSoln_[ gtspSoln_.size()-1 ]);
      }    
    }
    printf("\n");
  }

}

bool ODTSP::Dopt::isReplica( const int & node ){
  if ( node >= nt_*ne_ ){
    return true;
  }
  return false;
}


double ODTSP::Dopt::computeCost(){
  printf("computeCost()\n");
  printf("number of nodes : %d \n", nn_);
  printf("number of targets : %d \n", nt_);
  printf("number of clusters : %d \n", nc_);
  printf("number of replicas : %d \n", nr_);
  double cost = 0;
  if ( sol_ == LKH_NO_CLUSTERS || sol_ == GLKH_NO_CLUSTERS  ){

    // create a verbose flag to print this to terminal  
    for(int i = 0; i < atspSoln_.size() -1 ; i++){
      cost += (double)(C_[gtspSoln_[i]][gtspSoln_[i+1]]-offset_)/precision_;
      printf("C[%i][%i] = %3.1f \n", gtspSoln_[i], gtspSoln_[i+1] , (double)(C_[gtspSoln_[i]][gtspSoln_[i+1]]-offset_)/precision_ );
    }
    cost += (double) ( C_[atspSoln_[atspSoln_.size()-1]][atspSoln_[0]] ) / precision_;
    printf("C[%i][%i] = %3.1f \n", atspSoln_[atspSoln_.size()-1], atspSoln_[0] , (double)(C_[atspSoln_[atspSoln_.size()-1]][atspSoln_[0]]-offset_)/precision_ );
    
    
  }
  else if ( sol_ == GLKH_INTERSECTING_CLUSTERS ){
    for(int i = 0; i < gtspSoln_.size(); i++){
      int fromNode, toNode;
      if ( isReplica( gtspSoln_[i] ) ){
        int curCluster = whichCluster( gtspSoln_[i] );
        int numElems = list_->get_numElems( nt_ + curCluster );
        int exitNode = gtspSoln_[i] + ne_*(numElems-1);
        fromNode =  exitNode;
      }
      else {
        fromNode = gtspSoln_[i];
      }
      if ( i == gtspSoln_.size()-1 ){
        toNode = gtspSoln_[0];
      }
      else{ 
        toNode = gtspSoln_[i+1];
      }
      printf("C[%i][%i] = %3.1f \n", fromNode, toNode , (double)(C_[fromNode][toNode]-offset_)/precision_ );
      cost += (double)(C_[fromNode][toNode]-offset_)/precision_;    
    }
  }
  printf("Total = %3.1f\n", cost);
  // write code to compute GLK

  // write stats file
  std::string statFile = name_;
  statFile.append(".stats");
  FILE * pFile;
  pFile = fopen(statFile.c_str(),"a");
  fprintf(pFile,"%f ",cost);
  fclose(pFile);
  return(cost);
}


void ODTSP::Dopt::printCostMatrix(){
  for(int i = 0; i < nn_; i++){
    for(int j = 0; j < nn_; j++){
      if ( C_[i][j]== maxVal_ ){
        printf("x\t, ");
      }  
      else if ( offset_ == 0 ){
        if ( C_[i][j] > M_ ){
          printf("%3.0f+M\t, ", (double)(C_[i][j] - M_)/(double)precision_ );
        }
        else if ( C_[i][j] == -precision_ ){
          printf("?\t,  ");
        }
        else{
          printf("%3.0f\t, ", (double)(C_[i][j])/(double)precision_ );
        }
      } 
      else {
        if ( C_[i][j] > M_+offset_ ){
          printf("%3.0f+M+F\t, ", (double)(C_[i][j] - M_ - offset_)/(double)precision_ );
        }
        else if ( C_[i][j] == -precision_ ){
          printf("?\t,  ");
        }
        else if ( C_[i][j] >= offset_ ){
          printf("%3.0f+F\t, ", (double)((C_[i][j])-offset_)/(double)precision_ );
        }
        else{
          printf("%3.0f\t, ", (double)( C_[i][j] )/(double)precision_ );
        }
      }
      

    }
    printf("\n");
  }
}

void ODTSP::Dopt::printCostMatrixRaw(){
  for(int i = 0; i < nn_; i++){
    for(int j = 0; j < nn_; j++){
      printf("%3.0f\t, ", (double)((C_[i][j]))/(double)precision_ );
    }
    printf("\n");
  }
}


int ODTSP::Dopt::whichCluster( const int & node ){
  //printf("*** checking which cluster node %i is in ... \n",node);
  if ( node >= nt_*ne_ ){
    int pastNodes = nt_*ne_-1;
    for ( int i = 0; i < nc_; i++){
      pastNodes += clusterSize_[i]*ne_;
      //printf(" iter i=%i, pastNodes = %i \n",i,pastNodes);
      if ( node <= pastNodes ){
        return i;
      }      
    }
  }
  return -1;
}

bool ODTSP::Dopt::areReplicasSameCluster( const int & iNode, const int & jNode ){
  // both i and j must be >= nt_*ne_ to be replicas
  if ( iNode >= nt_*ne_ && jNode >= nt_*ne_ ){
    // check if they are part of the same cluster
    int iCluster = whichCluster( iNode );
    int jCluster = whichCluster( jNode );
    if ( iCluster == jCluster ){
      return true;
    }
  }
  return false;
}

void ODTSP::Dopt::printSoln(){

  // set up the file names
  std::string fName = name_;
  fName.append(".sol");
  // write the .sol file
  int from_list_indx,from_entryNode_indx,to_list_indx,to_entryNode_indx;
  double from_theta,to_theta,orbit;
  FILE *pFile;
  pFile = fopen(fName.c_str(),"w");
  resizeSolnVectors( gtspSoln_.size() );
  for(int i = 0; i < gtspSoln_.size(); ++i){
    if ( sol_ == GLKH_INTERSECTING_CLUSTERS ){
      from_list_indx = nodeListIndex_[ gtspSoln_[i] ];
      from_entryNode_indx = gtspSoln_[i]%ne_;
      if(i != gtspSoln_.size() - 1) {
        to_list_indx = nodeListIndex_[ gtspSoln_[i+1] ];
        to_entryNode_indx = gtspSoln_[i+1]%ne_;
      } 
      else {
        to_list_indx = nodeListIndex_[ gtspSoln_[0] ];
        to_entryNode_indx = gtspSoln_[0]%ne_;
      } 
    }
    else {
      from_list_indx = floor(gtspSoln_[i]/ne_);
      from_entryNode_indx = gtspSoln_[i]%ne_;
      if(i != gtspSoln_.size() - 1) {
        to_list_indx = floor(gtspSoln_[i+1]/ne_);
        to_entryNode_indx = gtspSoln_[i+1]%ne_;
      } 
      else {
        to_list_indx = floor(gtspSoln_[0]/ne_);
        to_entryNode_indx = gtspSoln_[0]%ne_;
      } 
    }
    
    // temp variables
    double fromX, fromY, fromHrad;
    double toX, toY, toHrad;
    double exitAngle, dir;

    // get the "from" point (the x,y coordiantes are NOT scaled up)
    dir = (double)orbitDir_[from_entryNode_indx];
    exitAngle = entryAngleRad_[from_entryNode_indx] 
                + dir*( list_->get_orbitAngleRad(from_list_indx) ) ;  
    transformEntryAngle( list_->get_orbitRadius(from_list_indx), exitAngle,dir, 
                         list_->get_centerX(from_list_indx), 
                         list_->get_centerY(from_list_indx),
                         fromX, fromY, fromHrad);
    // get the "to" point
    dir = (double)orbitDir_[to_entryNode_indx];
    exitAngle = entryAngleRad_[to_entryNode_indx] 
                + dir*( list_->get_orbitAngleRad(to_list_indx) ) ;         
    transformEntryAngle( list_->get_orbitRadius(to_list_indx), 
                         entryAngleRad_[to_entryNode_indx], dir, 
                         list_->get_centerX(to_list_indx), 
                         list_->get_centerY(to_list_indx),
                         toX, toY, toHrad);

    fprintf(pFile,"%d %d %f %f %f %f %f %f\n",from_list_indx,
            to_list_indx, fromX, fromY, fromHrad, toX, toY, toHrad);

    // save data
    fromX_[i] = fromX;
    fromY_[i] = fromY;
    fromHrad_[i] = fromHrad;
    fromCluster_[i] = from_list_indx;
    fromNode_[i] = from_entryNode_indx;
   

    toX_[i] = toX;
    toY_[i] = toY;
    toHrad_[i] = toHrad;
    toCluster_[i] = to_list_indx;
    toNode_[i] = to_entryNode_indx;
    
  }
  fprintf(pFile,"EOF");
  fclose(pFile);

}


ODTSP::Path ODTSP::Dopt::get_path(){
  list_->set_orbitPropertiesExit( fromCluster_, fromX_, fromY_, fromHrad_ );
  ODTSP::Path path(prob_, list_);
  path.set_numWpts(200);
  path.set_saveWaypoints(true);
  path.set_saveDubinsPaths(true);
  path.plan();  
  return path; 
}

void ODTSP::Dopt::print(){
  printf("--------------- DOPT SOLN -----------------------------\n");
  printf(" From                    | To           \t  | Cost \t | GTSP \n");
  printf(" No.  X,  Y,  H,         | No. X, Y, H  \t  |      \t | \n");
  int fromGtspNode;
  int toGtspNode;
  for (int i = 0; i < nl_; i++){
    fromGtspNode = gtspSoln_[i];
    if ( i == nl_-1 ){
      toGtspNode = gtspSoln_[0];
    } 
    else {
      toGtspNode = gtspSoln_[i+1];
    }
    printf("%i) %4.1f, %4.1f, %4.1f \t | %i) %4.1f, %4.1f, %4.1f, \t | %4.1f \t | %i->%i \n",fromCluster_[i],
    fromX_[i],fromY_[i],fromHrad_[i]*180.0/M_PI,
    toCluster_[i],toX_[i],toY_[i],toHrad_[i]*180.0/M_PI,
    ((double)C_[fromGtspNode][toGtspNode] - (double)M_)/ (double)precision_,
    fromGtspNode,toGtspNode );
  }
}

void ODTSP::Dopt::printGraphProperties(){
 printf("Number of nodes total : %i \n", nn_);
 printf("Number of entry states : %i \n", ne_);
 printf("Size of the list : %i \n", nl_);
 printf("Number of single targets : %i \n", nt_);
 printf("Number of clusters : %i \n", nc_);
 printf("Number of replicas : %i \n", nr_);
}


