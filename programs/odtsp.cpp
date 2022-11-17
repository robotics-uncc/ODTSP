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

#include<ODTSP_Solver.h>

bool plotPath( ODTSP::Path* path ){
  bool plotFlag = false;
  std::cout << " Plot? Requires Octave (sudo apt-get install octave) (0-no, 1-yes) : "; std::cin  >> plotFlag;
  // for plotting
  if ( plotFlag ){
    std::string fileName = "ODTSP.m";
    path->writePlotCommands(fileName);
    MathTools::runOctaveScript(fileName);
  }
  return true;
}

bool writePathXY( ODTSP::Path* path, std::string outputFile ){
  std::vector<double> x = path->get_xOverallPath();
  std::vector<double> y = path->get_yOverallPath();
  //std::vector<double> h = path->get_hOverallPath();
  FILE * fid = fopen( outputFile.c_str() , "a" );
  fprintf(fid,"pathXYH = [ ...\n");
  for (int i = 0; i < x.size(); i++ ){
    fprintf(fid,"%6.6f %6.6f;\n", x[i],y[i] ); //,h[i]); 
  }
  fprintf(fid,"];\n");
  fclose(fid);
}

int main(int argc, char** argv) {
  printf(" ******************************************************\n");
  printf("    ____     ______     ________    _____   _____      \n");
  printf("   / __ \\   (_  __ \\   (___  ___)  / ____\\ (  __ \\     \n");
  printf("  / /  \\ \\    ) ) \\ \\      ) )    ( (___    ) )_) )    \n");
  printf(" ( ()  () )  ( (   ) )    ( (      \\___ \\  (  ___/     \n");
  printf(" ( ()  () )   ) )  ) )     ) )         ) )  ) )        \n");
  printf("  \\ \\__/ /   / /__/ /     ( (      ___/ /  ( (         \n");
  printf("   \\____/   (______/      /__\\    /____/   /__\\        \n");
  printf("                                                       \n");
  printf(" Orbiting Dubins Traveling Salesman Problem           \n");
  printf(" A. Wolek, J. McMahon, Dec-2018                       \n");
  printf(" ******************************************************\n");
  printf(" Available Algorithms: \n");
  printf(" \t 0 - SYM_ANGLES_ETSP \n");
  printf(" \t 1 - SYM_ANGLES_ETSP_CURRENTS \n");
  printf(" \t 2 - NO_CLUSTER_GTSP \n");
  printf(" \t 3 - NO_CLUSTER_GTSP_CURRENTS \n");
  printf(" \t 4 - SIMPLE_CLUSTER_GTSP \n");
  printf(" \t 5 - SIMPLE_CLUSTER_GTSP_CURRENTS \n");
  printf(" \t 6 - EXHAUSTIVE_CLUSTER_GTSP \n");
  printf(" \t 7 - EXHAUSTIVE_CLUSTER_GTSP_CURRENTS \n");
  printf(" \t 8 - MAXIMAL_CLUSTER_TUNING_GTSP \n");
  printf(" \t 9 - MAXIMAL_CLUSTER_TUNING_GTSP_CURRENTS \n\n");
  printf(" \t 10 - MAXIMAL_CLUSTER_FEASIBLE_DUBINS \n\n");
  printf(" Available GTSP Solvers: \n");
  printf(" \t 0 - LKH \n");
  printf(" \t 1 - GLKH \n");
  printf("\nFor On-Screen Instructions: $ odtsp \n");
  printf("For Algorithms 0-1: $ odtsp file.scn alg file.sol \n");
  printf("For Algorithms 2-7: $ odtsp file.scn alg numEntryAngles tspSolver file.sol \n\n");
  printf("For Algorithms 8-10: $ odtsp file.scn alg numEntryAngles tspSolver tuningParam file.sol \n\n");
  printf("System Call: ");
  for ( int i = 0; i < argc; i++ ){
    printf("%s ",argv[i]);
  }
  printf("\n");
  // parameters
  int numTargs, precision;
  int tspSolverInt;
  double boxDim;
  ODTSP::RequiredParams params;
  double vehicleSpeed, windMagnitude, windDirRad, orbitAngle;
  double orbitRadiusMax, aspectAngleMaxRad, aspectAngleMaxDeg;
  std::vector<double> xTargs, yTargs;
  double windDirDeg;
  int alg = -1;
  int defaultScn = -1;
  int motionModelInt; 
  int numEntryAngles = -1;
  double tuningParam = 1;
  std::string sceneFilename;
  std::string outputFilename = "defaultOutput.sol";
  ODTSP::GTSP_SOLVER gtspSolver;

  // get scene
 if ( argc == 1 ){    
    std::cout << "Make a selection (0-9): "; std::cin  >> alg;    
    std::cout << "Use default scence (0-no,1-yes): "; std::cin  >> defaultScn;
    if ( defaultScn == 1 ){
      // set params      
      params.xTargs = {491.168, 55.546, 358.432, 333.343, 88.990, 177.341, 194.904, 292.540, 97.304, 86.564  };
      params.yTargs = {458.937, 177.273, 34.063, 98.208, 251.772, 196.381, 452.654, 121.693, 116.443, 8.033  };
      params.orbitAngleRad = 180.0 * M_PI /180.0;
      params.minTurnRadius = 30.0;
      params.orbitRadius = 60.0;
      orbitRadiusMax = 120;
      aspectAngleMaxRad = 30.0 * M_PI /180.0;
      vehicleSpeed = 1;
      windMagnitude = 0.5;
      windDirRad = 0;
    }
    else {
      std::cout << "Enter the dimensions of the box : "; std::cin  >> boxDim; 
      std::cout << "Enter the no. of targets : "; std::cin  >> numTargs;
      // set the target locations (a bounding box)
      double xmin = 0; double xmax = boxDim;
      double ymin = 0; double ymax = boxDim;
      if ( numTargs > 0 || boxDim <= 0 ){
        params.xTargs = MathTools::randomNumberUniformDistribution( xmin, xmax, numTargs); 
        params.yTargs = MathTools::randomNumberUniformDistribution( ymin, ymax, numTargs);
      }
      else{
        printf("Invalid entry.\n");
        return 0;
      }
      std::cout << "Enter the orbit angle (deg) : "; std::cin  >> orbitAngle; 
      params.orbitAngleRad = orbitAngle * M_PI /180.0;

      std::cout << "Enter the min. orbit radius (m) : "; std::cin  >> params.orbitRadius; 
      if ( alg >=4 && alg <= 9 ){
        std::cout << "Enter the max. orbit radius (m) : "; std::cin  >> orbitRadiusMax; 
        std::cout << "Enter the max. aspect angle (deg) : "; std::cin  >> aspectAngleMaxDeg;
        aspectAngleMaxRad = aspectAngleMaxDeg * M_PI / 180.0;
      }
      std::cout << "Enter motion model (0 - dubins, 1 - convected dubins ) : ";  
      std::cin  >> motionModelInt;
      if ( motionModelInt == 0){
        std::cout << "Enter the turn radius (m) : "; std::cin  >> params.minTurnRadius; 
      }
      else if ( motionModelInt == 1 ){
        std::cout << "Enter the turn radius (m) : "; std::cin  >> params.minTurnRadius; 
        std::cout << "Enter the vehicle speed (m/s) : "; std::cin  >> vehicleSpeed; 
        std::cout << "Enter the wind speed (m/s) : "; std::cin  >> windMagnitude; 
        std::cout << "Enter the wind direction (deg) : "; std::cin  >> windDirDeg; 
        windDirRad = windDirDeg * M_PI / 180.0;
      }
      else {
        printf("Invalid motion model. \n");
        return 0;
      }
      if ( alg > 1 ){ 
        std::cout << "Enter GTSP solver (0 - LKH, 1 - GLKH ) : ";  
        std::cin  >> tspSolverInt;
        if ( tspSolverInt == 0 ){
          gtspSolver = ODTSP::GTSP_SOLVER::LKH;
        }
        else if ( tspSolverInt == 1 ){
          gtspSolver = ODTSP::GTSP_SOLVER::GLKH;
        }
        else {
          printf("Invalid GTSP solver. \n");
          return 0;
        }
      }
      if ( alg > 7 ){
        std::cout << "Enter a orbit tuning parameter (0,1] : ";  
        std::cin  >> tuningParam;
        if ( tuningParam <= 0 || tuningParam > 1 ){
          printf("Invalid tuning param. \n");
          return 0;
        }
      }
    }

  }
  // argc = 4 , For Algorithms 0-1: $ odtsp file.scn alg file.sol 
  // argc = 6 , For Algorithms 2-7: $ odtsp file.scn alg numEntryAngles tspSolver file.sol 
  // argc = 7 , For Algorithms 8-9: $ odtsp file.scn alg numEntryAngles tspSolver tuningParam file.sol 
  else if ( argc == 4 || argc == 6 || argc == 7 ){ 
    printf("Received 4, 6, or 7 arguments\n");
    sceneFilename = argv[1];
    alg = atoi(argv[2]);    
    if ( argc == 4 ){
      outputFilename = argv[3];  
    }
    else if ( argc == 6 ){
      numEntryAngles = atoi(argv[3]);
      tspSolverInt = atoi(argv[4]);    
      outputFilename = argv[5];      
    }
    else if ( argc == 7 ){
      numEntryAngles = atoi(argv[3]);
      tspSolverInt = atoi(argv[4]); 
      tuningParam = atof(argv[5]);   
      outputFilename = argv[6];  
    }
    if ( tspSolverInt == 0 ){
      gtspSolver = ODTSP::GTSP_SOLVER::LKH;
    }
    else if ( tspSolverInt == 1 ){
      gtspSolver = ODTSP::GTSP_SOLVER::GLKH;
    }
    printf("\n");
    printf("Parsing %s...\n",sceneFilename.c_str());
    FILE * fid = fopen(sceneFilename.c_str(),"r");
    std::string curLine;
    std::ifstream sceneFile ( sceneFilename );

    while( getline(sceneFile, curLine) ) {
      if ( curLine.compare("ORBIT_ANGLE_DEG")==0 ){
        getline(sceneFile, curLine);
        params.orbitAngleRad = std::stod(curLine) * M_PI / 180.0;
        printf("ORBIT_ANGLE_DEG\n%3.3f\n",std::stod(curLine));
      }
      else if ( curLine.compare("ORBIT_RADIUS_MIN")==0 ){
        getline(sceneFile, curLine);
        params.orbitRadius = std::stod(curLine);
        printf("ORBIT_RADIUS_MIN\n%3.3f\n",params.orbitRadius);
      }
      else if ( curLine.compare("ORBIT_RADIUS_MAX")==0 ){
        getline(sceneFile, curLine);
        orbitRadiusMax = std::stod(curLine);
        printf("ORBIT_RADIUS_MAX\n%3.3f\n",orbitRadiusMax);
      }
      else if ( curLine.compare("POINTING_ERROR_MAX_DEG")==0 ){
        getline(sceneFile, curLine);
        aspectAngleMaxRad = std::stod(curLine) * M_PI / 180.0;
        printf("POINTING_ERROR_MAX_DEG\n%3.3f\n",std::stod(curLine));
      }
      else if ( curLine.compare("MIN_TURN_RADIUS")==0 ){
        getline(sceneFile, curLine);
        params.minTurnRadius = std::stod(curLine);
        printf("MIN_TURN_RADIUS\n%3.3f\n",params.minTurnRadius);
      }
      else if ( curLine.compare("VEHICLE_SPEED")==0 ){
        getline(sceneFile, curLine);
        vehicleSpeed = std::stod(curLine);
        printf("VEHICLE_SPEED\n%3.3f\n",vehicleSpeed);
      }
      else if ( curLine.compare("WIND_SPEED")==0 ){
        getline(sceneFile, curLine);
        windMagnitude = std::stod(curLine);
        printf("WIND_SPEED\n%3.3f\n",windMagnitude);
      }
      else if ( curLine.compare("WIND_DIR_DEG")==0 ){
        getline(sceneFile, curLine);
        windDirRad = std::stod(curLine) * M_PI / 180.0;
        printf("WIND_DIR_DEG\n%3.3f\n", std::stod(curLine) );
      }
      else if ( curLine.compare("PRECISION")==0 ){
        getline(sceneFile, curLine);
        precision = std::stoi(curLine);
        printf("PRECISION\n%d\n",precision);
      }
      else if ( curLine.compare("BOX_DIM")==0 ){
        getline(sceneFile, curLine);
        boxDim = std::stod(curLine);
        printf("BOX_DIM\n%3.3f\n",boxDim); 
      }
      else if ( curLine.compare("NUM_TARGS")==0 ){
        getline(sceneFile, curLine);
        numTargs = std::stoi(curLine);
        printf("NUM_TARGS\n%d\n",numTargs);
      }
      else if ( curLine.compare("TARGS")==0 ){
        getline(sceneFile, curLine);
        printf("TARGS\n");
        float x, y;
        while( curLine.compare("EOF")!=0 ){
          sscanf(curLine.c_str(),"%f%f",&x, &y);
          xTargs.push_back( x );
          yTargs.push_back( y );          
          printf("%3.4f %3.4f\n", xTargs.back(), yTargs.back());
          getline(sceneFile, curLine);
        }
        printf("EOF\n");
      }
    }    
    fclose(fid);
    params.xTargs = xTargs;
    params.yTargs = yTargs;
  }
  else {
    printf("\nFor On-Screen Instructions: $ odtsp \n");
    printf("For Algorithms 0-1: $ odtsp file.scn alg file.sol \n");
    printf("For Algorithms 2-7: $ odtsp file.scn alg numEntryAngles file.sol \n");
    printf("For Algorithms 8-10: $ odtsp file.scn alg numEntryAngles tspSolver tuningParam file.sol \n\n");
  }
  ODTSP::Path* path;
  if ( alg < 0 || alg > 10 ){
    return 0;
  }
  if ( numEntryAngles == -1 && alg > 1 ){
    std::cout << "Number of Entry Angles? (2-16) : "; std::cin  >> numEntryAngles;
  }
  if ( alg == (int)ODTSP::Algorithm::SYM_ANGLES_ETSP ){
    // solve
    ODTSP::Solver solver( ODTSP::Algorithm::SYM_ANGLES_ETSP , params);
    path = solver.getPath();
    solver.writeSolnFile( outputFilename );
    writePathXY( path, outputFilename ); 
  }
  else if ( alg == (int)ODTSP::Algorithm::SYM_ANGLES_ETSP_CURRENTS ){
    // solve
    ODTSP::Solver solver( ODTSP::Algorithm::SYM_ANGLES_ETSP_CURRENTS , params, 
                          vehicleSpeed, windMagnitude, windDirRad );
    path = solver.getPath();
    solver.writeSolnFile( outputFilename );
    writePathXY( path, outputFilename ); 
  }
  else if ( alg == (int)ODTSP::Algorithm::NO_CLUSTER_GTSP ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, params, numEntryAngles );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );    
    writePathXY( path, outputFilename );  
  }
  else if ( alg == (int)ODTSP::Algorithm::NO_CLUSTER_GTSP_CURRENTS ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, params, numEntryAngles , 
                          vehicleSpeed, windMagnitude, windDirRad );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );    
    writePathXY( path, outputFilename ); 
  }
  else if ( alg == (int)ODTSP::Algorithm::SIMPLE_CLUSTER_GTSP || 
            alg == (int)ODTSP::Algorithm::EXHAUSTIVE_CLUSTER_GTSP ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, 
                          params, numEntryAngles, 
                          orbitRadiusMax, aspectAngleMaxRad );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );    
    writePathXY( path, outputFilename ); 
  }
  else if ( alg == (int)ODTSP::Algorithm::SIMPLE_CLUSTER_GTSP_CURRENTS || 
            alg == (int)ODTSP::Algorithm::EXHAUSTIVE_CLUSTER_GTSP_CURRENTS ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, 
                          params, numEntryAngles, 
                          orbitRadiusMax, aspectAngleMaxRad, 
                          vehicleSpeed, windMagnitude, windDirRad );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );  
    writePathXY( path, outputFilename ); 
    
  }
  else if ( alg == (int)ODTSP::Algorithm::MAXIMAL_CLUSTER_TUNING_GTSP ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, 
                          params, numEntryAngles, 
                          orbitRadiusMax, aspectAngleMaxRad, tuningParam );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );  
    writePathXY( path, outputFilename );     
  }
  else if ( alg == (int)ODTSP::Algorithm::MAXIMAL_CLUSTER_TUNING_GTSP_CURRENTS || 
            alg == (int)ODTSP::Algorithm::MAXIMAL_CLUSTER_FEASIBLE_DUBINS_GTSP ){
    ODTSP::Solver solver( (ODTSP::Algorithm)alg , gtspSolver, 
                          params, numEntryAngles, 
                          orbitRadiusMax, aspectAngleMaxRad, tuningParam ,
                          vehicleSpeed, windMagnitude, windDirRad );
    path = solver.getPath(); 
    solver.writeSolnFile( outputFilename );  
    writePathXY( path, outputFilename );     
  }
   
  return 0;
}
