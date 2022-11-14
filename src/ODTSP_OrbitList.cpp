#include "ODTSP_OrbitList.h"
#include<Miniball.hpp>
#include<algorithm>

// -----------------------------------------------------------------------------
// Initialization 
// -----------------------------------------------------------------------------

bool ODTSP::OrbitList::initialize( ODTSP::Problem* prob, 
                                   ODTSP::ClusterMode mode ){
  if ( prob->isDefined() ){
    prob_ = prob;
    mode_ = mode;
    generateTargetList();
    return true;
  }
  else {
    printf("ODTSP::OrbitList constructed with an ill-defined problem.\n");
    return false;
  }  
}

ODTSP::OrbitList::OrbitList( ODTSP::Problem* prob, 
                             ODTSP::ClusterMode mode ){
  initialize(prob, mode);
}


// -----------------------------------------------------------------------------
// High-level functionality 
// -----------------------------------------------------------------------------

void ODTSP::OrbitList::generateTargetList(){
  int n = prob_->get_numTargs();
  int numEdges;    
  vvi F(n);
  list_.resize(n);
  // intialize all of the targets before clustering
  for ( int i = 0; i < n; i++){
    vi elems;
    elems.push_back(i);
    // initialize the target list as all single types from the problem 
    (list_[i]).initialize( elems, prob_->get_xTargs(i), prob_->get_yTargs(i), 
                                                  prob_->get_orbitRadius(), 
                                                  prob_->get_orbitAngleRad(),
                                                  ODTSP::TargetType::Single );   
  }
  // check which clustering type is needed 
  if ( mode_ == ODTSP::ClusterMode::SimpleCluster ){
    getAdjacencyMatrix(F, numEdges, 4.0 );
    if ( numEdges > 0 ){      
      vvi comboSet;
      getNeighborSets(comboSet, F);
      deleteRepeats(comboSet);
      vvvi binnedCombos;
      binCombosBySize(binnedCombos, comboSet);
      breakAndAddUniqueClustersToList(binnedCombos);
    }
  }
  else if ( mode_ == ODTSP::ClusterMode::Exhaustive ){
    getAdjacencyMatrix(F, numEdges, 4.0 );
    vvi comboSet;   
    if ( numEdges > 0 ){
      // graph file populated by writeDIMACSfile function 
      std::string graphfile = "cliquerFile.txt";
      writeDIMACSfile( graphfile, F, n, numEdges );
      vvi comboSet = allCliques( graphfile );  
      runMiniballAndAddValidToList(comboSet);
    }
  }
  else if ( mode_ == ODTSP::ClusterMode::MaximalClusterTuning ){
    getAdjacencyMatrix(F, numEdges, tuningParam_ );
    if ( numEdges > 0 ){
      // graph file populated by writeDIMACSfile function 
      std::string graphfile = "cliquerFile.txt";
      writeDIMACSfile( graphfile, F, n, numEdges );
      vvi comboSet = maximalCliques( graphfile );
      vvvi binnedCombos;
      binCombosBySize(binnedCombos, comboSet);
      breakAndAddUniqueClustersToList(binnedCombos);
    }
    else{
      //printf("no clusters formed.\n");
    }
  }
}

// -----------------------------------------------------------------------------
// Simple Cluster 
// -----------------------------------------------------------------------------

// functions below are used for ODTSP::ClusterMode::SimpleCluster
void ODTSP::OrbitList::getNeighborSets(vvi & comboSet, vvi F){  
  // determine how many candidates are connected to each node 
  int n = prob_->get_numTargs();  // number of points  
  // for each node add connectivity vector to the comboSet
  // this allows us to remove repeats easily before binning them 
  for (int i = 0; i < n; i++){
    // create a connectivity vector based on the F matrix row from above 
    // (i.e., a vector with all potential candidates)
    vi curNodeConnectivity;   
    for (int j = 0; j < n; j++){
      if ( F[i][j] == 1 ){
        curNodeConnectivity.push_back(j); 
      }
    } 
    if ( curNodeConnectivity.size() > 1){
      comboSet.push_back( curNodeConnectivity );
    } 
  }  
}

// -----------------------------------------------------------------------------
// Exhaustive Cluster 
// -----------------------------------------------------------------------------

// functions below are used for ODTSP::ClusterMode::Exhaustive
void ODTSP::OrbitList::powerSetsOfMaximalCliques(vvi & comboSet){
  // comboset already contains all maximal cliques but they may overlap and
  // here we want all powersets
  for ( int i = 0; i < comboSet.size(); i++ ){
    vi combo = comboSet[i];
    if ( combo.size() > 2 ){
      for (int k = 2; k < combo.size()-1; k++){
        enumerateCombinations(combo, k, comboSet);
        deleteRepeats(comboSet);
      }    
      
    }
  }

}

// -----------------------------------------------------------------------------
// Shared functions 
// -----------------------------------------------------------------------------

void ODTSP::OrbitList::binCombosBySize(vvvi &binnedCombos, vvi &comboSet ){
  // determine how many elements are in each combo 
  int largestBin = 1;
  for ( int i = 0; i < comboSet.size(); i++){
    if ( (comboSet[i]).size() > largestBin ){
      largestBin = (int)(comboSet[i]).size();   
    }
  }

  // create a bin which stores cluster candidates of different sizes
  // the first entry binnedCombos[0] stores candidates of size numCandsMax
  // the second entry binnedCombos[1] stores candidates of size numCandsMax-1
  // ...
  // the last entry binnedCombos[numCandsMax-2] stores candidates of size 2 
  if ( largestBin > 1 ){
    binnedCombos.resize(largestBin - 1);
  }

  // bin results 
  for (int i = 0; i < comboSet.size(); i++){
    int sizeOfSet = (comboSet[i]).size();
    (binnedCombos[ largestBin - sizeOfSet ]).push_back( comboSet[i] );
  } 
}

void ODTSP::OrbitList::getAdjacencyMatrix(vvi &F, int & numEdges, double eta){
  int numVertices = prob_->get_numTargs();  // number of points
  double dsqrmax = eta*eta*prob_->get_dmax()*prob_->get_dmax();   
  for (int i = 0; i < numVertices; i++){
    (F[i]).resize(numVertices);
  } 
  numEdges = 0;
  // compute connectivity matrix if F[i][j] = 1 then it is possible
  // for node i to be clustered with node j
  for (int i = 0; i < numVertices; i++){ 
    //for (int j = 0; j < i; j++){
    for (int j = 0; j < numVertices; j++){
      // compute distance squared between pt. i and pt. j
      double dsqr = (prob_->get_xTargs(i) - prob_->get_xTargs(j))*
                    (prob_->get_xTargs(i) - prob_->get_xTargs(j)) + 
                    (prob_->get_yTargs(i) - prob_->get_yTargs(j))*
                    (prob_->get_yTargs(i) - prob_->get_yTargs(j)) ; 
      if ( dsqr <= dsqrmax ){
        F[i][j] = 1;
        if ( i != j ){
          numEdges = numEdges + 1;
        }
      }
      else {
        F[i][j] = 0;
      }
    }
  }
  numEdges = numEdges/2;
  return;
}

void ODTSP::OrbitList::breakAndAddUniqueClustersToList( vvvi & binnedCombos ){  
  int numBins = binnedCombos.size();
  vi usedVertices; // use this vector to store list of vertices that have 
                   // already been assigned to a valid cluster 
  // go through each bin 
  for ( int i = 0; i < numBins; i++ ){
    // go through each cluster within the bin 
    for (int j = 0; j < (binnedCombos[i]).size(); j++ ){
      vi candCluster = binnedCombos[i][j]; // cluster currently being considered
      // compute the difference between candCluster and usedValues
      vi uniqueCluster; // the resulting set from the difference 
      if ( !std::is_sorted( candCluster.begin(), candCluster.end() )){
        throw std::runtime_error("candCluster not sorted."); 
      }
      if ( !std::is_sorted( usedVertices.begin(), usedVertices.end() )){
        throw std::runtime_error("usedVertices not sorted."); 
      }  
      std::set_difference( candCluster.begin(), candCluster.end(),
                           usedVertices.begin(), usedVertices.end(),
                           std::back_inserter(uniqueCluster) );
      // if the unique cluster is not same size as the candidate, then place 
      // it into a lower bin to be processed later 
      if ( uniqueCluster.size() < candCluster.size() 
           && uniqueCluster.size() >= 2){
        (binnedCombos[numBins - uniqueCluster.size()]).push_back(uniqueCluster);
        //(binnedCombos[numBins - uniqueCluster.size() - 1]).push_back(uniqueCluster);
      }
      else {
        // check if this is a valid cluster 
        ODTSP::Target targ;
        if ( isValid(uniqueCluster, targ) && uniqueCluster.size() >= 2 ){
          list_.push_back( targ ); // final result goes here 
          usedVertices.insert( usedVertices.end(), uniqueCluster.begin(), 
                               uniqueCluster.end() ); 
          std::sort( usedVertices.begin(), usedVertices.end() ); 
        }     
        // if it is not valid, then remove a "distant" point and demote to a 
        // lower bin 
        else {
          removeDistantPoint( uniqueCluster );
          if ( uniqueCluster.size() < candCluster.size() 
               && uniqueCluster.size() >= 2){
            (binnedCombos[numBins - uniqueCluster.size()]).push_back(uniqueCluster);
            //(binnedCombos[numBins - uniqueCluster.size() - 1]).push_back(uniqueCluster);
          }  
        }
      }
    }
  }
  // now remove all of the ODTSP::TargetType::Single types that are in 
  // clusters, working backwards, delete the entries of list_
  // (this is possible because all the clusters have been 'push_back' to end)
  for (int i = usedVertices.size()-1; i >= 0; i = i - 1){
    list_.erase( list_.begin() + usedVertices[i] );
  }
}

void ODTSP::enumerateCombinations(const vi & baseCombo, 
                                  int k, 
                                  vvi & comboSet){
  int n = baseCombo.size();
  if ( k <= 1 ){
    throw std::runtime_error("ODTSP::OrbitList::enumerateCombinations: k <= 1\n"); 
  }
  else if (k > n - 1) {
    throw std::runtime_error("ODTSP::OrbitList::enumerateCombinations: k > n"); 
  }

  std::vector<bool> v(n); // a "selector" vector
  std::fill(v.begin(), v.begin() + k, true); // make the first k elements true
  bool permuteCycleIncomplete = true; 
  // keep adding combos until all permutations of selector vector considered 
  while ( permuteCycleIncomplete ){ 
    vi currentComb;  
    for (int i = 0; i < n; ++i) {
       if ( v.at(i) ) {
          currentComb.push_back( baseCombo.at(i) );
       }
    }
    std::sort(currentComb.begin(), currentComb.end()); // sort before appending 
    if ( currentComb.size() > 1 ){
      bool isRepeat = false;
      for ( int i = 0; i < comboSet.size(); i++ ){
        if ( equal_seq(currentComb, comboSet[i]) ){
          isRepeat = true;
        }
      }
      if ( !isRepeat ){
        comboSet.push_back(currentComb); // append       
      }
      //comboSet.shrink_to_fit();
    }
    // permute elements of v and check if permutation has completed full cycle
    permuteCycleIncomplete = std::prev_permutation(v.begin(), v.end()); 
  }
}

void ODTSP::OrbitList::deleteRepeats(vvi & comboSet){
  // we sort first before checking for uniqueness
  if ( comboSet.size() > 0 ){
    vi emptyElems;
    for (int i = 0; i < comboSet.size(); i++ ){
      //assert( (comboSet[i]).size() >= 1 );
      if ( (int)(comboSet[i]).size() < 1 ){
        emptyElems.push_back(i);
      }
    }
//    for (int i = 0; i < emptyElems.size(); i++ ){
//      int indexToErase = emptyElems[ emptyElems.size() - i - 1 ];
//      comboSet.erase( comboSet.begin() + indexToErase );
//    }
//    
    for (int i = 0; i < comboSet.size(); i++ ){
      assert( (comboSet[i]).size() >= 1 );
    }        
    assert( comboSet.size() >= 1 );
    std::sort(comboSet.begin(), comboSet.end(), ODTSP::compare_seq);
    vvi::iterator it;
    it = std::unique(comboSet.begin(), comboSet.end(), ODTSP::equal_seq);
    //comboSet.erase( it, comboSet.end() );
    comboSet.resize( std::distance(comboSet.begin(),it) ); 
    comboSet.shrink_to_fit();
  }
}

bool ODTSP::OrbitList::isValid( vi & candCluster, ODTSP::Target & targ ){
  if ( candCluster.size() > 1){
    // shorthand notation:
    typedef double* const* PointIterator; 
    typedef const double* CoordIterator;
    typedef Miniball::
      Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;
    // for each combo run the miniball algorithm to compute smallest-circle
    // center and radius
    int d = 2;  // dimension of position vector 
    int numPts = candCluster.size(); // number of points in cluster
    double** ap = new double*[numPts]; // matrix of points
    for (int j = 0; j < numPts; j++) {
      double* p = new double[d];
      p[0] = prob_->get_xTargs( candCluster[j] );
      p[1] = prob_->get_yTargs( candCluster[j] );
      ap[j] = p;
    }
    MB mb (d, ap, ap+numPts); // main run of miniball
    const double* center = mb.center();
    double clusterRadius = std::sqrt(mb.squared_radius() );
    // clean-up
    for (int k = 0; k < numPts; k ++){
      delete [] ap[k];
    }
    delete [] ap;
    // check if the cluster is still feasible with the miniball radius
    if ( clusterRadius <= prob_->get_dmax() ){
      // use the miniball result to specify the cluster center and radius
      targ.initialize( candCluster, center[0],     
                                    center[1], 
                                    clusterRadius + prob_->get_orbitRadius(), 
                                    prob_->get_orbitAngleRad(),
                                    ODTSP::TargetType::Cluster );
      return true;
    }
    return false;
  }
  else {
    return false;
  }
}

void ODTSP::OrbitList::removeDistantPoint( vi & cluster ){
  // only possible to remove points if there are at least 2
  if ( cluster.size() >= 2){
    // compute center of mass
    double xc = 0;
    double yc = 0;
    int numCands = cluster.size();
    for ( int i = 0; i < numCands; i++ ){
      xc = xc + prob_->get_xTargs( cluster[i] );
      yc = yc + prob_->get_yTargs( cluster[i] );
    }
    xc = xc/(double)numCands;
    yc = yc/(double)numCands;
    // compute distance from center of mass to each point in the clusterData
    vi d(numCands);  
    int dMaxInd = 0;
    double dMax = -1;
    for ( int i = 0; i < numCands; i++){
      double dsqr =  ( prob_->get_xTargs( cluster[i] )-xc )
                    *( prob_->get_xTargs( cluster[i] )-xc )  + 
                     ( prob_->get_yTargs( cluster[i] )-yc )
                    *( prob_->get_yTargs( cluster[i] )-yc ) ;
      if ( dsqr > dMax ){
        dMaxInd = i;  
        dMax = dsqr;
      }
    }
    // remove the point with maximum distance from the cluster
    cluster.erase( cluster.begin() + dMaxInd ); 
  }
}


// -----------------------------------------------------------------------------
// Sequence utils  
// -----------------------------------------------------------------------------

// Is 'second' greater than 'first' ?
bool ODTSP::compare_seq(vi & first, vi & second){
//  printf("comparing sequence first:\n");
//  printClusterVector(first);
//  printf("comparing sequence second:\n");
//  printClusterVector(second);

  int numDigitsFirst = first.size();
  int numDigitsSecond = second.size();
 // printf(" (numDigitsFirst, numDigitsSecond) = (%d, %d)\n", numDigitsFirst, numDigitsSecond);
  int minDigits = std::min(numDigitsFirst, numDigitsSecond);
  if ( minDigits == 0 ){
    if (numDigitsFirst > numDigitsSecond){ 
      return false;
    }
    else {
      return true;
    }
  }
  else{
    for (int i = 0; i < minDigits; i++){
      if ( second[i] > first[i] ){
        return true;
      }
      else if ( second[i] < first[i] ){
        return false;
      }
    }
    if (numDigitsFirst > numDigitsSecond){ 
      return false;
    }
    else {
      return true;
    }
  }
}

bool ODTSP::OrbitList::validSequence(vi seq){
  // check that the sequence provided is the correct size
  if ( seq.size() > (list_).size() ){
    throw std::runtime_error("OrbitList::validSequence: sequence size greater than list."); 
  }
  // sort a copy of the input sequence
  vi seqSorted = seq;
  std::sort( seqSorted.begin(), seqSorted.end() );
  // check that the elements of the input sequence are within range of the list
  if ( seqSorted[0] < 0 || seqSorted[seqSorted.size()-1] > (list_).size()-1 ){
    throw std::runtime_error("OrbitList::validSequence: sequence elements outside of list range."); 
  }
  // check that the input sequence is unique
  auto last = std::unique( seqSorted.begin(), seqSorted.end() );
  //seqSorted.resize( std::distance(seqSorted.begin(),last) ); 
  if ( last!=seqSorted.end() ){
    throw std::runtime_error("OrbitList::validSequence: sequence elements not unique."); 
  }
}

bool ODTSP::OrbitList::sortAndTruncate(vi seq){
  // check if sequence is valid
  validSequence( seq );
  // create the sorted list by extracting elements from list_ and make index vec
  std::vector<ODTSP::Target> listSorted;
  for (int i = 0; i < seq.size(); i++){
    listSorted.push_back( list_[ seq[i] ] );
  }
  // redefine list_
  list_ = listSorted;
  return true;
}

bool ODTSP::OrbitList::elemsInCommon(int i, int j, vi & elems){ 
    if ( elems.size() > 0){
      // the elems vector should be empty when passed to this function, 
      // otherwise it will always return true
      throw std::runtime_error("ODTSP::OrbitList::elemsInCommon: input vi elems already defined.");
    }
    vi v1 = (list_[i]).get_elems();
    vi v2 = (list_[j]).get_elems();
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(elems));  
    if ( elems.size() > 0){
      return true;
    }
    else{
      return false;
    }
}

bool ODTSP::OrbitList::elemsInCommon(int i, int j){ 
  vi elems;
  return elemsInCommon(i,j,elems);
}


bool ODTSP::equal_seq (std::vector<int> & first, std::vector<int> & second){
  if ( first.size() == second.size() ){
    for (int i = 0; i < first.size(); i++){
      if ( first.at(i) != second.at(i) ){
        return false;
      }
    }
  }
  else {
    return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
// Miniball interface 
// -----------------------------------------------------------------------------

void ODTSP::OrbitList::runMiniballAndAddValidToList(const vvi & comboSet){
  if ( comboSet.size() > 0 ){
    // shorthand notation:
    typedef double* const* PointIterator; 
    typedef const double* CoordIterator;
    typedef Miniball::
      Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;
    // for each combo run the miniball algorithm to compute smallest-circle
    // center and radius    
    for ( int i = 0; i < comboSet.size(); i++ ){
      int d = 2;  // dimension of position vector 
      int numPts = (comboSet[i]).size(); // number of points in cluster
      double** ap = new double*[numPts]; // matrix of points
      for (int j = 0; j < numPts; j++) {
        double* p = new double[d];
        p[0] = prob_->get_xTargs( (comboSet[i])[j] );
        p[1] = prob_->get_yTargs( (comboSet[i])[j] );
        ap[j] = p;
      }
      MB mb (d, ap, ap+numPts); // main run of miniball
      const double* center = mb.center();
      double clusterRadius = std::sqrt(mb.squared_radius() );
      // check if the cluster is still feasible with the miniball radius
      if ( clusterRadius <= prob_->get_dmax() ){     
        ODTSP::Target targ; 
        // use the miniball result to specify the cluster center and radius
        targ.initialize( comboSet[i], center[0],     
                                      center[1], 
                                      clusterRadius + prob_->get_orbitRadius(), 
                                      prob_->get_orbitAngleRad(),
                                      ODTSP::TargetType::Cluster );
        list_.push_back(targ); // add to list 
      }
      else{
        //printf("miniball says cluster not valid. Cluster radius : %3.3f, but dmax = %3.3f \n", clusterRadius, prob_->get_dmax() );
      }
      // clean-up
      for (int k = 0; k < numPts; k ++){
        delete [] ap[k];
      }
      delete [] ap;
    }
  }
  else{
    printf("Warning: comboSet.size() = 0\n");
  }
}


// -----------------------------------------------------------------------------
// Cliquer interface 
// -----------------------------------------------------------------------------

void ODTSP::OrbitList::writeDIMACSfile( std::string graphfile, vvi &F, int numVert, int numEdge ){
  FILE *pFile;
  pFile = fopen(graphfile.c_str(),"w");
  // Line: p FORMAT NODES EDGES
  // Each line contains one 'p' line which describes the dimension of the graph
  // FORMAT is used for consistency with older formats and should contain the word 'edge'
  // NODES is the number of vertices (numbered from 1 to n))
  // EDGES is the number of edges 
  // Cliquer ignores the FORMAT and EDGES entries, but they must be there anyway
  fprintf(pFile,"p edge %d %d\n", numVert, numEdge );
  // go through each vertex 
  for (int i = 0; i < numVert; i++){
    // since the matrix is symemtric we only go through upper right triangle
    for (int j = 0; j < i; j++){
      // if an edge is preset, add an 'e' line to the file
      if ( F[i][j] == 1 && i != j){
        // Line: e W V 
        // Specifies that there is an edge between vertices W and V. This line is *not* repeated
        // as "e V W"
        fprintf(pFile,"e %d %d\n",i+1,j+1);
      }
    }
  }
  fclose(pFile);
}

// request only the maximal cliques for the graph 
vvi ODTSP::OrbitList::allCliques( std::string graphfile ){
  // Warning! This requires that the "cl" cliquer executable location is in your PATH variable 
  std::string command = "cl --all -m 2 ";
  std::string outputFile = graphfile;
  outputFile.append(".out");
  command.append( graphfile.c_str() );
  command.append(" > ");
  command.append( outputFile.c_str() );
  // send the command to Cliquer
  system(command.c_str());
  // open the output
  FILE *pFile;
  pFile = fopen(outputFile.c_str(),"r");
  // exmaple output:
  // size=4, weight=4:   1 3 5 8
  int cliqueSize, cliqueWeight;
  std::string cliqueString;
  vvi cliqueSet;
  int node;
  std::string curLine;
  std::ifstream sceneFile ( outputFile );
  std::string delimeter = " ";
  while( getline(sceneFile, curLine) ) {  
    sscanf(curLine.c_str(),"size=%d, weight=%d:   ", &cliqueSize, &cliqueWeight);
    //printf("size=%d, weight=%d:   \n", cliqueSize, cliqueWeight);
    vi clique(cliqueSize); 
    // find the colon position
    std::string msg = curLine;
    std::size_t valStart = 0;
    std::size_t valEnd = msg.find(":"); // find first occurance of delimeter
    valStart = valEnd + 3;
    valEnd = msg.find(delimeter,valStart);
    int i = 0;
    while ( valEnd != std::string::npos ){ // while there is a value to extract
      // get the substring for the current number being parsed
      std::string curValString = msg.substr(valStart,valEnd - valStart);
      int val; 
      if ( !curValString.empty() && curValString.compare(delimeter)!=0 ){       
        //printf("curValString : %s \n", curValString.c_str() );
        val = std::stoi(curValString); // convert string to double
        clique[i] = val-1; // add value to the result
        i = i + 1;
      }
      // check if there are other remaining space delimited pairs
      valStart = valEnd+1;
      valEnd = msg.find(delimeter,valEnd+1);
    }
    if ( valStart < valEnd){ // there is still a final number to be extracted
      // get the substring for the current number being parsed
      std::string curValString = msg.substr(valStart,valEnd - valStart);
      int val; 
      // if string is not empty
      if ( !curValString.empty() && curValString.compare(delimeter)!=0 ){       
        val = std::stoi(curValString); // convert string to double
        clique[i] = val-1;
      }
      // add value to the result    
    }
    //    printf("found clique: ");
    //    for ( int j = 0; j < clique.size(); j++ ){
    //      printf(" %d, ", clique[j]);
    //    }
    //    printf("\n");
    // sort (needed for comparisons later)
    std::sort( clique.begin(), clique.end() );
    cliqueSet.push_back(clique);   
  }
  fclose(pFile);
  return cliqueSet;
}


// request only the maximal cliques for the graph 
vvi ODTSP::OrbitList::maximalCliques( std::string graphfile ){
  // Warning! This requires that the "cl" cliquer executable location is in your PATH variable 
  std::string command = "cl --all -m 2 -x ";
  std::string outputFile = graphfile;
  outputFile.append(".out");
  command.append( graphfile.c_str() );
  command.append(" > ");
  command.append( outputFile.c_str() );
  // send the command to Cliquer
  system(command.c_str());
  // open the output
  FILE *pFile;
  pFile = fopen(outputFile.c_str(),"r");
  // exmaple output:
  // size=4, weight=4:   1 3 5 8
  int cliqueSize, cliqueWeight;
  std::string cliqueString;
  vvi cliqueSet;
  int node;
  std::string curLine;
  std::ifstream sceneFile ( outputFile );
  std::string delimeter = " ";
  while( getline(sceneFile, curLine) ) {  
    sscanf(curLine.c_str(),"size=%d, weight=%d:   ", &cliqueSize, &cliqueWeight);
    //printf("size=%d, weight=%d:   \n", cliqueSize, cliqueWeight);
    vi clique(cliqueSize); 
    // find the colon position
    std::string msg = curLine;
    std::size_t valStart = 0;
    std::size_t valEnd = msg.find(":"); // find first occurance of delimeter
    valStart = valEnd + 3;
    valEnd = msg.find(delimeter,valStart);
    int i = 0;
    while ( valEnd != std::string::npos ){ // while there is a value to extract
      // get the substring for the current number being parsed
      std::string curValString = msg.substr(valStart,valEnd - valStart);
      int val; 
      if ( !curValString.empty() && curValString.compare(delimeter)!=0 ){       
        //printf("curValString : %s \n", curValString.c_str() );
        val = std::stoi(curValString); // convert string to double
        clique[i] = val-1; // add value to the result
        i = i + 1;
      }
      // check if there are other remaining space delimited pairs
      valStart = valEnd+1;
      valEnd = msg.find(delimeter,valEnd+1);
    }
    if ( valStart < valEnd){ // there is still a final number to be extracted
      // get the substring for the current number being parsed
      std::string curValString = msg.substr(valStart,valEnd - valStart);
      int val; 
      // if string is not empty
      if ( !curValString.empty() && curValString.compare(delimeter)!=0 ){       
        val = std::stoi(curValString); // convert string to double
        clique[i] = val-1;
      }
      // add value to the result    
    }
    // sort (needed for comparisons later)
    std::sort( clique.begin(), clique.end() );
    cliqueSet.push_back(clique);   
  }
  fclose(pFile);
  return cliqueSet;
}


// -----------------------------------------------------------------------------
// set/get functions 
// -----------------------------------------------------------------------------

bool ODTSP::OrbitList::set_orbitProperties( vi & seq, 
                                            vd & entryAngles, 
                                            vi & orbitDirs ){
  // check if sequence is valid
  validSequence( seq );
  // set the entry angle and orbit dir of each target
  for ( int i = 0; i < seq.size(); i++){
    (list_[seq[i]]).set_entryAngleOrbitDir( entryAngles[i], orbitDirs[i] );
  }
  sortAndTruncate(seq);
  return true;
}

bool ODTSP::OrbitList::set_orbitProperties( vi & seq, vd & entryX, 
                                            vd & entryY, 
                                            vd & entryHrad ){
  // check if sequence is valid
  validSequence( seq );
  // set the entry state of each target
  for ( int i = 0; i < seq.size(); i++){  
    (list_[seq[i]]).set_entryState( entryX[i], entryY[i], entryHrad[i]);
  }
  sortAndTruncate(seq);
  return true;
}

bool ODTSP::OrbitList::set_orbitPropertiesExit( vi & seq, 
                                                vd & exitX, 
                                                vd & exitY, 
                                                vd & exitHrad ){
  // check if sequence is valid
  validSequence( seq );
  // set the exit state of each target
  for ( int i = 0; i < seq.size(); i++){
    printf(" i : %i of %i \n", i , seq.size() );
    (list_[seq[i]]).print();
    (list_[seq[i]]).set_exitState( exitX[i], exitY[i],  exitHrad[i]);
  }
  sortAndTruncate(seq);
  return true;
}

bool ODTSP::OrbitList::validIndex(int ind){
  if ( ind >= 0 && ind <= (int)list_.size() ){
    return true;
  }
  printf(" requested ind = %i, but list size is %i \n",ind, (int)list_.size() );
  throw std::runtime_error("ODTSP::OrbitList::validIndex: requsted index invalid.");
  return false;
}

double ODTSP::OrbitList::get_exitX(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_exitX();
  }
}

double ODTSP::OrbitList::get_exitY(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_exitY();
  }
}

double ODTSP::OrbitList::get_exitHrad(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_exitHrad();
  }
}

double ODTSP::OrbitList::get_exitAngleRad(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_exitAngleRad();
  }
}


double ODTSP::OrbitList::get_entryX(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_entryX();
  }
}

double ODTSP::OrbitList::get_entryY(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_entryY();
  }
}

double ODTSP::OrbitList::get_entryHrad(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_entryHrad();
  }
}


double ODTSP::OrbitList::get_entryAngleRad(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_entryAngleRad();
  }
}

int ODTSP::OrbitList::get_orbitDir(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_orbitDir();
  }
}

double ODTSP::OrbitList::get_orbitRadius(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_orbitRadius();
  }
}

double ODTSP::OrbitList::get_orbitAngleRad(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_orbitAngleRad();
  }
}

double ODTSP::OrbitList::get_centerX(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_centerX();
  }
}

double ODTSP::OrbitList::get_centerY(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_centerY();
  }
}

ODTSP::TargetType ODTSP::OrbitList::get_type(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_type();
  }
}

int ODTSP::OrbitList::get_numElems(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_numElems();
  }
}

vi ODTSP::OrbitList::get_elems(int ind){
  if ( validIndex(ind) ){ 
    return (list_[ind]).get_elems();
  }
}


bool ODTSP::OrbitList::contains_elem( int listInd, 
                                      int elemInd ){
  return (list_[listInd]).contains_elem(elemInd);
}

bool ODTSP::OrbitList::contains_elem( int listInd, 
                                      int elemInd,
                                      int & entryInd ){
  return (list_[listInd]).contains_elem(elemInd, entryInd);
}

// -----------------------------------------------------------------------------
// Print/Display 
// -----------------------------------------------------------------------------

void ODTSP::printClusterVector( const vi & cluster ){
  for (int j = 0; j < cluster.size(); j++){
    printf("%i , ", cluster[j]);
  }
  printf("\n");
}


void ODTSP::printComboSet( const vvi & comboSet ){
  for (int i = 0; i < comboSet.size(); i++){
    printf("comboSet[%i] : ",i);
    printClusterVector( comboSet[i] );
  }
}

void ODTSP::printBinnedCombos( const vvvi & binnedCombos ){
  for ( int i = 0; i < binnedCombos.size(); i++ ){
    printf("** bin %i **\n",i);
    printComboSet( binnedCombos[i] );
  }
}

void ODTSP::OrbitList::print(){
  for (int i = 0; i < list_.size(); i++){
    printf("--- Orbit List ID: %d --- \n",i);
    (list_[i]).print();
  }
  if ( list_.size() == 0 ){
    printf("Orbit list_size() = 0 \n");
  }
}

void ODTSP::OrbitList::printClusters(){
  for (int i = 0; i < list_.size(); i++){
    if ( (list_[i]).get_type() == ODTSP::TargetType::Cluster ){
      printf("--- Orbit List ID: %d --- \n",i);
      (list_[i]).print();
    }
  }
}
