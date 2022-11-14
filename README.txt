 ******************************************************
    ____     ______     ________    _____   _____      
   / __ \   (_  __ \   (___  ___)  / ____\ (  __ \     
  / /  \ \    ) ) \ \      ) )    ( (___    ) )_) )    
 ( ()  () )  ( (   ) )    ( (      \___ \  (  ___/     
 ( ()  () )   ) )  ) )     ) )         ) )  ) )        
  \ \__/ /   / /__/ /     ( (      ___/ /  ( (         
   \____/   (______/      /__\    /____/   /__\        
                                                       
 Orbiting Dubins Traveling Salesman Problem           
 A. Wolek, J. McMahon, Sep-2018                       
 ******************************************************
 Available Algorithms: 
 	 0 - SYM_ANGLES_ETSP 
 	 1 - SYM_ANGLES_ETSP_CURRENTS 
 	 2 - NO_CLUSTER_GTSP 
 	 3 - NO_CLUSTER_GTSP_CURRENTS 
 	 4 - SIMPLE_CLUSTER_GTSP 
 	 5 - SIMPLE_CLUSTER_GTSP_CURRENTS 
 	 6 - EXHAUSTIVE_CLUSTER_GTSP 
 	 7 - EXHAUSTIVE_CLUSTER_GTSP_CURRENTS 

----------------------------------------------------------

To install dependencies, and build this project for the first time:

Move into your desired directory:
$ cd <path-to-your-desired-directory>

Download the project:
git clone https://github.com/arturwolek/optimalRID.git

First install/build libODTSP:
$ cd optimalRID
$ cd libODTSP
$ . install_deps.sh 

This should complete the install and build. 

----------------------------------------------------------

To modify and re-build project:

Move into the libODTSP folder
$ cd libODTSP

Clean all the build files
$ . clean.sh -all

Build
$ . build.sh


